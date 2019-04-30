
import sys
import argparse
import uuid
import subprocess



def shell_cmd(
	cmd, 
	calling_name):

	'''
	Wrapper for subprocess.check_output().

	Attempts calling shell command. Prints
	error if any errors are raised and exits
	program if so.

	Parameters
	----------

	cmd : str
		String containing full shell command
	calling_name : 
		Name of function (or other environment) 
		from which shell_cmd was called - 
		useful for debugging.
	'''

	try:
		subprocess.check_output(
        	cmd, 
        	shell = True, 
        	stderr = subprocess.STDOUT)

	except subprocess.CalledProcessError as e:
		print e.output
		sys.exit("Subprocess command failed in " + 
        	     calling_name + ". Exiting . . . ")




def call_bedtools_sort(
	input_file,
	outdir,
	outfile = None,
	bedtools_exec = "bedtools",
	remove_original = True):

	'''
	Sort input bed file or bg file
	using bedtools sort with the following
	command:

	$ bedtools sort -i input_file > outdir/outfile_name

	Parameters
	----------

	input_file : str
		Full path to input file
	outdir : str
		Full path to output directory
	outfile : str 
		If set, provides full path for sorted file.
		If not set, function will insert ".sorted"
		prior to the input_file suffix.
	bedtools_exec : str
		Path to bedtools binary. Default
		"bedtools" assumes bedtools in user's
		PATH.
	remove_original : bool
		Whether to remove the original (unsorted file).
		Often there is no reason to keep it, and as such
		default is True.

	'''

	if not outfile:

		input_file_split = input_file.split(".")
		input_file_prefix = ".".join(input_file_split[0:-1])
		suffix = input_file_split[-1]

		outfile = (input_file_prefix + 
				   ".sorted." + 
				   suffix)

	bedtools_sort_command = " ".join([
		bedtools_exec,
		"sort",
		"-i",
		input_file,
		">",
		outfile])

	shell_cmd(bedtools_sort_command,
			  "call_bedtools_sort - sort")

	if remove_original:

		rm_command = " ".join([
			"rm",
			input_file])

		shell_cmd(rm_command,
			      "call_bedtools_sort - rm")

	return outfile


def expand_bedgraph(
	path_to_bg,
	outfile = None,
	remove_original = True):

	'''
	Expands a bedgraph file such that intervals with
	length > 1 become multiple intervals with
	length = 1. Writes the expanded file to disk.
	This expansion facilitates using single-base resolution. 

	Parameters
	----------
	path_to_bg : str
		Path to bedgraph file
	outfile : 
		Path to output file.  If unset the original file
		path with "expanded" inserted prior to the suffix is
		used.
	remove_original : bool
		Whether to remove the original (unsorted file).
		Often there is no reason to keep it, and as such
		default is True.
	'''


	if not outfile:

		input_file_split = path_to_bg.split(".")
		input_file_prefix = ".".join(input_file_split[0:-1])
		suffix = input_file_split[-1]

		outfile = (input_file_prefix + 
				   ".expanded." + 
				   suffix)



	expanded_file = open(outfile, 'w')

	with open(path_to_bg, 'r') as file:

		for line in file:

			entry = line.strip().split()

			start = int(entry[1])
			end = int(entry[2])

			gap = end - start

			if gap == 1:

				expanded_file.write(line)

			else:

				for i in range(start, end):

					expanded_file.write(
						"\t".join([entry[0], 
								   str(i), 
								   str(i+1)] + entry[3:]) + "\n")

	if remove_original:

		rm_command = " ".join([
			"rm",
			path_to_bg])

		shell_cmd(rm_command,
			      "expand_bedgraph - rm")

	expanded_file.close()


	return outfile



def call_bedtools_intersect(
	path_to_region_bed,
	path_to_bg,
	outdir,
	filename_prefix = "",
	bedtools_exec = "bedtools"):

	'''
	Calls bedtools intersect using subprocess
	module with the following effective command:

	$ bedtools intersect -a path_to_region_bed \
						 -b path_to_bg \
						 -wa \
						 -wb \
						 -sorted \
						 > \
						 output_file

	The idea is to retrieve bedgraph information
	(e.g. phastCons scores, mapability, etc) for
	features of interest (e.g. regions surrounding
	splice sites).  Note that input bed AND begraph
	(bg) files MUST be sorted - it is easiest to sort 
	them using sortBed/bedtools sort.  Failure to do
	so will lead to unexpected behavior.


	Parameters
	----------

	path_to_region_bed : str
		Path to region bed file
	path_to_bg : str
		Path to bedgraph file
	outdir : str
		Path to output directory
	filename_prefix : str
		Any prefix to be added to filename - e.g.
		"human" for human output if process will 
		be repeated for multiple species
	bedtools_exec : str
		Path to bedtools binary. Default
		"bedtools" assumes bedtools in user's
		PATH.
	'''

	if filename_prefix != "":

		filename_prefix = "_" + filename_prefix

	initial_intersect_outpath = (outdir + 
								 "/" + 
								 str(uuid.uuid4()) + 
								 filename_prefix + 
								 "_bg_overlap.bedlike")

	bedtools_intersect_command = " ".join([
		bedtools_exec,
		"intersect",
		"-a",
		path_to_region_bed,
		"-b",
		path_to_bg,
		"-wa",
		"-wb",
		"-sorted",
		">",
		initial_intersect_outpath])

	shell_cmd(bedtools_intersect_command,
			  "call_bedtools_intersect")

	return initial_intersect_outpath 




def process_overlap(
	initial_intersect_outpath,
	outdir,
	name_field_header,
	region_length,
	filename_prefix = "",
	remove_original = True):

	'''
	Writes simplified output file based on 
	bedtools intersect.

	Assumes adjacent regions appear sequentially

	Parameters
	----------

	initial_intersect_outpath : str
		Path to bedtools intersect output file
	outfile_path : str
		Path to output file
	name_field_header : list
		list if strings constituting the header
		for information contained in field 3 of
		the bedtools intersect output
	prefix : str
		Prefix to add to output filename
	region_length : int
		Size of the regions
	remove_original : bool
		Whether to remove the original (unsorted file).
		Often there is no reason to keep it, and as such
		default is True.
	'''

	if filename_prefix != "":

		filename_prefix = filename_prefix + "_"


	outfile_path = outdir + "/" + filename_prefix + "metaregion_scores.tsv"

	output_file = open(outfile_path, 'w')

	header = "\t".join([
		"shared_coordinate",
		"score"
		] + name_field_header)

	output_file.write(header + "\n")

	with open(initial_intersect_outpath, 'r') as file:

		full_identifier = ""

		for line in file:
			entry = line.split()

			strand = entry[5]
			score = entry[9]

			identifier_split = entry[3].split(";")
			host_feature_id = identifier_split[-1]

			if entry[3] != full_identifier:
				full_identifier = entry[3]

				if strand == "+":
					counter = 0

				elif strand == "-":
					counter = region_length - 1

			output_entry = "\t".join([
				str(counter),
				score] + 
				identifier_split)

			output_file.write(output_entry + "\n")

			if strand == "+":

				counter += 1

			elif strand == "-":

				counter -= 1


	if remove_original:

		rm_command = " ".join([
			"rm",
			initial_intersect_outpath])

		shell_cmd(rm_command,
			      "process_overlap - rm")

	output_file.close()



def main():

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--region_bed",
		type = str,
		help = "Path to region bed file",
		required = True)

	parser.add_argument(
		"--bedgraph",
		type = str,
		help = "Path to begGraph file",
		required = True)

	parser.add_argument(
		"--outdir",
		type = str,
		help = "Path to output directory",
		required = True)

	parser.add_argument(
		"--name_fields",
		type = str,
		help = "Comma separated list of names of name field elements",
		required = True)

	parser.add_argument(
		"--region_length",
		type = int,
		help = "Number of bases in each region",
		required = True)

	parser.add_argument(
		"--sort_bed",
		action = "store_true",
		help = "Sort bed file if not sorted")

	parser.add_argument(
		"--sort_bedgraph",
		action = "store_true",
		help = ("Sort bedgraph file if not sorted. " + 
			   "If used, this will likely take a while. " + 
			   "Set if not already sorted"))

	parser.add_argument(
		"--expand_bedgraph",
		action = "store_true",
		help = ("Expand bedgraph to 1-line per base - " + 
				" set this if not already expanded"))

	parser.add_argument(
		"--bedtools_exec",
		type = str,
		help = "Path to bedtools binary. default = 'bedtools'",
		default = "bedtools")

	parser.add_argument(
		"--prefix",
		type = str,
		help = "Prefix for output_files, default = ''",
		default = "")


	args = parser.parse_args()


	if args.sort_bed:

		region_bed = call_bedtools_sort(
						args.region_bed,
						args.outdir,
						outfile = None,
						bedtools_exec = args.bedtools_exec,
						remove_original = True)

	else:

		region_bed = args.region_bed



	if args.sort_bedgraph:

		bedgraph = call_bedtools_sort(
						args.bedgraph,
						args.outdir,
						outfile = None,
						bedtools_exec = args.bedtools_exec,
						remove_original = True)

	else:

		bedgraph = args.bedgraph


	if args.expand_bedgraph:

		expanded_bedgraph = expand_bedgraph(
						bedgraph,
						outfile = None,
						remove_original = True)

	else:

		expanded_bedgraph = bedgraph


	initial_intersect_outpath = call_bedtools_intersect(
									region_bed,
									expanded_bedgraph,
									args.outdir,
									filename_prefix = args.prefix,
									bedtools_exec = args.bedtools_exec)


	process_overlap(
		initial_intersect_outpath,
		args.outdir,
		args.name_fields.split(","),
		args.region_length,
		filename_prefix = args.prefix,
		remove_original = True)




if __name__ == "__main__":

	main()


				
						

			









	

