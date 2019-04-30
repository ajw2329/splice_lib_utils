#!/usr/bin/python

import sys
from splice_lib import splice_lib
import copy
import argparse
import subprocess
from maxentpy import maxent # must be installed
import gen_methods
import pandas as pd



def get_donors_acceptors(event):
	'''
	Identify donor, acceptor splice sites in a 
	splice_lib-style event dictionary entry

	Parameters
	----------
	event : dict
	    dictionary entry from splice_lib 
	    standard event dict (or similar)
	    containing alternative splicing
	    event details

	Returns
	-------
	included_donors : dict
		contains dict with keys (int) corresponding
		to included form donor coords, and 
		values (int) corresponding to matching 
		acceptor coords
	included_acceptors : dict
		same as above but with keys (int) corresponding
		to acceptors and values (int) corresponding to
		matching donors
	excluded_donors : dict
		same as included_donors but for excluded form
	excluded_acceptors: dict
		same as included_acceptors but for excluded form
	donor_index : int
		index of donor coord in junction tuple
		0 if strand == "+", 1 if strand == "-"
	acceptor_index :  int
		index of acceptor coord in junction tuple
		1 if strand == "+", 0 if strand == "-"	
	rev : bool	
		Logical indicating whether splice site order
		should be reversed (convenient for later 
		functions when strand == "-"). True if
		strand == "-", False if strand == "+".
	'''

	strand = event["strand"]

	if strand == "+":

		donor_index = 0
		acceptor_index = 1
		rev = False

	elif strand == "-":

		donor_index = 1
		acceptor_index = 0
		rev = True

	else:

		sys.exit("Unexpected strand character in event entry", strand, event)

	included_donors = {}
	included_acceptors = {}

	excluded_donors = {}
	excluded_acceptors = {}

	for junction in event["included_junctions"]:

		included_donors[int(junction[donor_index])] = int(junction[acceptor_index])
		included_acceptors[int(junction[acceptor_index])] = int(junction[donor_index])

	for junction in event["excluded_junctions"]:

		excluded_donors[int(junction[donor_index])] = int(junction[acceptor_index])
		excluded_acceptors[int(junction[acceptor_index])] = int(junction[donor_index])

	return (included_donors, 
		    included_acceptors, 
		    excluded_donors, 
		    excluded_acceptors, 
		    donor_index, 
		    acceptor_index, 
		    rev)



def get_competitive_junctions(included_donors, 
	                          included_acceptors, 
	                          excluded_donors, 
	                          excluded_acceptors):
	'''
	Identify competitive junctions in splicing 
	event, based on dicts producted by calling
	get_donors_acceptors(event).  

	Competitive junctions are those in which either
	the donor or the acceptor is present in both the
	included and excluded isoforms of the alternative
	splicing event. If that is true, then the 
	isoform-specific splice sites are potentially
	"competing" for access to the pan-isoform splice site.
	Note that this categorization may fail in contexts 
	where the splicing decision may be forced by 
	prior decisions (e.g. in the case of alt-TSS where
	the alternative donor may in one isoform not be
	transcribed at all).

	Parameters
	----------
	included_donors : dict
		contains dict with keys (int) corresponding
		to included form donor coords, and 
		values (int) corresponding to matching 
		acceptor coords
	included_acceptors : dict
		same as above but with keys (int) corresponding
		to acceptors and values (int) corresponding to
		matching donors
	excluded_donors : dict
		same as included_donors but for excluded form
	excluded_acceptors: dict
		same as included_acceptors but for excluded form


	Returns
	-------
	included_competitive_junctions : list
		list of tuples each having length 2
		in which the first element is the 
		donor coordinate (int) of the competitive
		junction and the second element is the acceptor 
		coordinate.
	excluded_competitive_junctions : list
		lists of tuples each having length 2
		in which the first element is the 
		donor coordinate (int) of the competitive
		junction and the second element is the acceptor 
		coordinate.

	'''

	included_competitive_junctions = []
	excluded_competitive_junctions = []

	for donor, acceptor in included_donors.items():

		if (donor in excluded_donors) or (acceptor in excluded_acceptors):

			included_competitive_junctions.append((donor, acceptor))

	for donor, acceptor in excluded_donors.items():

		if (donor in included_donors) or (acceptor in included_acceptors):

			excluded_competitive_junctions.append((donor, acceptor))

	return (included_competitive_junctions, 
		    excluded_competitive_junctions)




def get_competing_splice_sites(event, 
	                           acceptor_coordinates, 
	                           donor_coordinates):
	'''
	Modifies event dict entry in place (no returns).
	Identifies splice sites unique to specific
	alternative splicing event isoforms (i.e.
	splice sites that belong to competitive junctions
	and are unique to a specific isoform). Adds
	donor and acceptor dicts to event dict entry
	(keys "donor" and "acceptor" respectively), which
	each include values describing coordinates, the 
	isoform to which the splice site is unique,
	coordinates of the partner (e.g. acceptor if the
	discussing the "donor" dict), and the coordinates
	of the competitor (i.e. the splice site in the
	other isoform that is spliced to the same partner),
	and the index of the competitive junction in the 
	isoform when they are all ordered from 5' -> 3'.
	Acceptor and donor coordinates are also added to
	sets acceptor_coordinates and donor_coordinates
	as described below.

	Parameters
	----------
	event : dict
	    dictionary entry from splice_lib 
	    standard event dict (or similar)
	    containing alternative splicing
	    event details
	acceptor_coordinates : set
		empty set which will be populated by
		strings describing individual acceptor
		coordinates in the following format:
		chromosome + "_" + coordinate + "_" + strand
		e.g. "chr1_12573532_+"
	donor_coordinates : set
		empty set which will be populated by
		strings describing individual donor
		coordinates in the following format:
		chromosome + "_" + coordinate + "_" + strand
		e.g. "chr2_31907124_-"


	Modifications to input
	----------------------

	event : dict
		dict is modified as described above,
		to contain new values
		"competitive_unique_donors" : dict
			dict with coords of competitive donors
			unique to a specific isoform as keys,
			and dicts as values each with the 
			following values : 

				"acceptor" : int
					coordinate of the acceptor to
					which the donor is spliced
				"unique_form" : str
					event isoform to which the donor
					is unique ("included" or "excluded")
				"competing_donor" : int
					coordinate of the donor in the other
					isoform with which the key donor
					competes
				"competitive_junction_index" : int
					0-based index of the junction
					(amongst the other competitive
					junctions of the same isoform 
					sorted from 5' -> 3')

	acceptor_coordinates : set
		acceptor coordinates are added to this set
		as strings containing chromosome and strand
		e.g. "chr1_12573532_+"
	donor_coordinates : set
		donor coordinates are added to this set
		as strings containing chromosome and strand
		e.g. "chr1_12573532_+"

		
	'''

	strand = event["strand"]
	chrom = event["chrom"]

	(included_donors, 
     included_acceptors, 
     excluded_donors, 
     excluded_acceptors, 
     donor_index, 
     acceptor_index, 
     rev) = get_donors_acceptors(event)

	(included_competitive_junctions, 
	 excluded_competitive_junctions) = get_competitive_junctions(included_donors, 
	                                                             included_acceptors, 
	                                                             excluded_donors, 
	                                                             excluded_acceptors)

	donors_in_included_competitive_junctions = []
	acceptors_in_included_competitive_junctions = []

	donors_in_excluded_competitive_junctions = []
	acceptors_in_excluded_competitive_junctions = []


	competitive_unique_donors = {}
	competitive_unique_acceptors = {}

	for i, junction in enumerate(sorted(included_competitive_junctions, 
										reverse = rev)):

		donor = int(junction[0])
		acceptor = int(junction[1])

		donors_in_included_competitive_junctions.append(donor)
		acceptors_in_included_competitive_junctions.append(acceptor)

		if donor not in excluded_donors:

			try:
				competitive_unique_donors[donor] = {"acceptor": acceptor,
												    "unique_form": "included",
												    "competing_donor": excluded_acceptors[acceptor],
												    "competitive_junction_index": i + 1
												    }
			except KeyError as e:

				print "KeyError of interest!"
				print e
				print included_competitive_junctions
				print included_donors
				print included_acceptors
				print excluded_donors
				print excluded_acceptors
				print event
				sys.exit()



		if acceptor not in excluded_acceptors:

			competitive_unique_acceptors[acceptor] = {"donor": donor,
												      "unique_form": "included",
												      "competing_acceptor": excluded_donors[donor],
												      "competitive_junction_index": i + 1
												      }




	for i, junction in enumerate(excluded_competitive_junctions):

		donor = int(junction[0])
		acceptor = int(junction[1])

		donors_in_excluded_competitive_junctions.append(donor)
		acceptors_in_excluded_competitive_junctions.append(acceptor)

		if donor not in included_donors:

			competitive_unique_donors[donor] = {"acceptor": acceptor,
											    "unique_form": "excluded",
											    "competing_donor": included_acceptors[acceptor],
											    "competitive_junction_index": i + 1
											    }

		if acceptor not in included_acceptors:

			competitive_unique_acceptors[acceptor] = {"donor": donor,
											          "unique_form": "excluded",
											          "competing_acceptor": included_donors[donor],
											          "competitive_junction_index": i + 1
											          }


	event["competitive_unique_donors"] = copy.deepcopy(competitive_unique_donors)
	event["competitive_unique_acceptors"] = copy.deepcopy(competitive_unique_acceptors)


	for donor in set(donors_in_included_competitive_junctions + 
					 donors_in_excluded_competitive_junctions):

		basic_coords = chrom + "_" + str(donor) + "_" + strand

		donor_coordinates.add(basic_coords)



	for acceptor in set(acceptors_in_included_competitive_junctions + 
						acceptors_in_excluded_competitive_junctions):

		basic_coords = chrom + "_" + str(acceptor) + "_" + strand

		acceptor_coordinates.add(basic_coords)





def write_splice_site_region_bed(
	event_dict,
	outdir,
	species,
	exonic = 30,
	intronic = 150):

	bedfile_out = open(
		outdir + "/" + species + "_splice_site_regions.bed", 
		'w')

	for event, event_val in event_dict.iteritems():

		donor_coords = event_val.get("competitive_unique_donors")
		acceptor_coords = event_val.get("competitive_unique_acceptors")

		all_coords = {}
		strand = event_val["strand"]
		event_type = event_val["event_type"]
		chrom = event_val["chrom"]		

		if donor_coords:

			for donor, donor_val in donor_coords.iteritems():

				all_coords[donor] = "donor"
				all_coords[donor_val["competing_donor"]] = "donor"		
				all_coords[donor_val["acceptor"]] = "acceptor"

		if acceptor_coords:

			for acceptor, acceptor_val in acceptor_coords.iteritems():

				all_coords[acceptor] = "acceptor"
				all_coords[acceptor_val["competing_acceptor"]] = "acceptor"		
				all_coords[acceptor_val["donor"]] = "donor"

		if strand == "+":

			rev = False

		elif strand == "-":

			rev = True

		sorted_coords = sorted(all_coords.keys(), reverse = rev)

		for index, i in enumerate(sorted_coords):

			if ((all_coords[i] == "donor" and strand == "+") or 
				(all_coords[i] == "acceptor" and strand == "-")):

				start = i - exonic
				end = i + intronic

			else:

				start = i - intronic
				end = i + exonic

			entry = [chrom, 
					 str(start), 
					 str(end), 
					 event_type + ";" + str(index) + ";" + event,
					 "1000",
					 strand]

			bedfile_out.write("\t".join(entry) + "\n")

	bedfile_out.close()




def write_splice_site_bed(splice_sites, 
						  splice_site_type, 
						  outdir, 
						  species, 
						  donor_exonic = 3, 
						  donor_intronic = 6, 
						  acceptor_exonic = 3, 
						  acceptor_intronic = 20):
	'''
	
	Parameters
	----------

	splice_sites : set
		set containing splice site coords as strings
		in the following format: "chr1_23123123_+"
	splice_site_type : str
		either "donor" or "acceptor"
	outdir : str
		path to output directory
	donor_exonic : int
		Number of exonic bases to include in donor coords
	donor_intronic : int
		Number of intronic basis to include in donor coords
	acceptor_exonic : int
		Number of exonic bases to include in donor coords
	acceptor_intronic : int
		Number of intronic basis to include in donor coords

	Returns
	-------

	None

	Instead, outputs a bed file containing coordinates
	'''

	bedfile = open(outdir + "/" + species + "_" + splice_site_type + ".bed", 'w')

	for splice_site in splice_sites:

		splice_site_split = splice_site.split("_")


		chrom = splice_site_split[0]
		coord = int(splice_site_split[1])
		strand = splice_site_split[2]

		if strand == "+":

			if splice_site_type == "donor":

				start = coord - donor_exonic + 1 - 1
				end = coord + donor_intronic

			elif splice_site_type == "acceptor":

				start = coord - acceptor_intronic - 1
				end = coord + donor_exonic - 1

		elif strand == "-":

			if splice_site_type == "donor":

				start = coord - donor_intronic - 1
				end = coord + donor_exonic - 1

			elif splice_site_type == "acceptor":

				start  = coord - acceptor_exonic + 1 - 1
				end = coord + acceptor_intronic

		bedfile.write("\t".join([chrom, 
								 str(start), 
								 str(end), 
								 splice_site, 
								 "1000",
								 strand]) + "\n")

	bedfile.close()

	#for minus strand: end = donor + 2, start = donor - 6 (for 1-based)
	# iiiiii[E]EE
	# for plus strand: start = donor - 2, end = donor + 6
	# EE[E]iiiiii
	# 12 3 456789

	#for plus strand: end = donor + 2, start = donor - 20 (for 1-based)
	# iiiiiiiiiiiiiiiiiiii[E]EE
	# for minus strand: start = donor - 2, end = donor + 20
	# EE[E]iiiiiiiiiiiiiiiiiiii



def getfasta(outdir, 
			 genome_fasta, 
			 species, 
			 bedtools_path):

	donor_cmd = [bedtools_path,
				 "getfasta", 
				 "-bed", 
				 outdir + "/" + species + "_donor.bed",
				 "-tab",
				 "-fi",
				 genome_fasta,
				 "-fo",
				 outdir + "/" + species + "_donor.fastatab",
				 "-s",
				 "-name"]

	acceptor_cmd = [bedtools_path,
					"getfasta", 
					"-bed", 
					outdir + "/" + species + "_acceptor.bed",
					"-tab",
					"-fi",
					genome_fasta,
					"-fo",
					outdir + "/" + species + "_acceptor.fastatab",
					"-s",
					"-name"]

	gen_methods.run_bash_cmd(donor_cmd)
	gen_methods.run_bash_cmd(acceptor_cmd)


def read_and_score_fasta(outdir, 
	                     species, 
	                     donor_dinucleotide_start = 3, 
	                     acceptor_dinucleotide_start = 18):

	donor_dict = {}
	acceptor_dict = {}

	acceptor_scorefile = open(outdir + "/" + species + "_acceptor_scores.tsv", 'w')
	acceptor_scorefile.write("\t".join(["splice_site_type", 
		                                "location",
		                                "seq",
		                                "score",
		                                "dinucleotide",
		                                "dinucleotide_is_standard"]) + "\n")

	donor_scorefile = open(outdir + "/" + species + "_donor_scores.tsv", 'w')
	donor_scorefile.write("\t".join(["splice_site_type", 
		                             "location",
		                             "seq",
		                             "score",
		                             "dinucleotide",
		                             "dinucleotide_is_standard"]) + "\n")

	with open(outdir + "/" + species + "_donor.fastatab", 'r') as file:

		donor_matrix = maxent.load_matrix5()

		for line in file:

			entry = line.strip().split("\t")

			key = entry[0].split("(")[0]
			seq = entry[1].upper()
			dinucleotide = seq[donor_dinucleotide_start:donor_dinucleotide_start+2]
			standard_dinucleotide = dinucleotide == "GT"


			donor_dict[key] = {"seq": seq, 
			                   "score": maxent.score5(seq, donor_matrix) if "N" not in seq else "NA", 
			                   "dinucleotide": dinucleotide, 
			                   "standard_dinucleotide": standard_dinucleotide}

			donor_scorefile.write("\t".join(["donor", 
											 key, 
											 seq, 
											 str(donor_dict[key]["score"]), 
											 dinucleotide, 
											 str(standard_dinucleotide)]) + "\n")

			

	with open(outdir + "/" + species + "_acceptor.fastatab", 'r') as file:

		acceptor_matrix = maxent.load_matrix3()

		for line in file:

			entry = line.strip().split("\t")

			key = entry[0].split("(")[0]
			seq = entry[1].upper()
			dinucleotide = seq[acceptor_dinucleotide_start:acceptor_dinucleotide_start+2]
			standard_dinucleotide = dinucleotide == "AG"

			acceptor_dict[key] = {"seq": seq, 
			                      "score": maxent.score3(seq, acceptor_matrix) if "N" not in seq else "NA", 
			                      "dinucleotide": dinucleotide, 
			                      "standard_dinucleotide": standard_dinucleotide}

			acceptor_scorefile.write("\t".join(["acceptor", 
				                                key, 
				                                seq, 
				                                str(acceptor_dict[key]["score"]), 
				                                dinucleotide, 
				                                str(standard_dinucleotide)]) + "\n")

	donor_scorefile.close()
	acceptor_scorefile.close()

	return donor_dict, acceptor_dict
			

def associate_seqs_scores_with_events(event_dict, 
	                                  donor_dict, 
	                                  acceptor_dict, 
	                                  outdir, 
	                                  species,
	                                  full_command_line):
	'''
		event is the dictionary entry associated with an event in event_dict
	'''

	events = []
	splice_site_type = []
	unique_form = []
	competitive_junction_index = []
	primary_key = []
	primary_seq = []
	primary_score = []
	primary_dinucleotide = []
	primary_dinucleotide_standard = []
	competitor_key = []
	competitor_seq = []
	competitor_score = []
	competitor_dinucleotide = []
	competitor_dinucleotide_standard = []
	partner_key = []
	partner_seq = []
	partner_score = []
	partner_dinucleotide = []
	partner_dinucleotide_standard = []

	for event, event_val in event_dict.iteritems():

		chrom = event_val["chrom"]
		strand = event_val["strand"]

		if "competitive_unique_acceptors" in event_val:

			for acceptor, acceptor_val in event_val["competitive_unique_acceptors"].iteritems():

				events += [event]

				competitive_junction_index += [acceptor_val["competitive_junction_index"]]

				acceptor_key = chrom + "_" + str(acceptor) + "_" + strand
				competing_acceptor_key = chrom + "_" + str(acceptor_val["competing_acceptor"]) + "_" + strand
				donor_key = chrom + "_" + str(acceptor_val["donor"]) + "_" + strand

				splice_site_type += ["acceptor"]
				primary_key += [acceptor_key]
				competitor_key += [competing_acceptor_key]
				partner_key += [donor_key]
				unique_form += [acceptor_val["unique_form"]]


				if acceptor_key in acceptor_dict:

					acceptor_val["acceptor_score"] = acceptor_dict[acceptor_key]["score"]
					acceptor_val["acceptor_seq"] = acceptor_dict[acceptor_key]["seq"]
					acceptor_val["acceptor_key"] = acceptor_key
					acceptor_val["acceptor_dinucleotide"] = acceptor_dict[acceptor_key]["dinucleotide"]
					acceptor_val["acceptor_standard_dinucleotide"] = acceptor_dict[acceptor_key]["standard_dinucleotide"]

					primary_seq += [acceptor_dict[acceptor_key]["seq"]]
					primary_score += [acceptor_dict[acceptor_key]["score"]]
					primary_dinucleotide += [acceptor_dict[acceptor_key]["dinucleotide"]]
					primary_dinucleotide_standard += [acceptor_dict[acceptor_key]["standard_dinucleotide"]]

				else:

					primary_seq += ["NA"]
					primary_score += ["NA"]
					primary_dinucleotide += ["NA"]
					primary_dinucleotide_standard += ["NA"]					

				if competing_acceptor_key in acceptor_dict:

					acceptor_val["competing_acceptor_score"] = acceptor_dict[competing_acceptor_key]["score"]
					acceptor_val["competing_acceptor_seq"] = acceptor_dict[competing_acceptor_key]["seq"]
					acceptor_val["competing_acceptor_key"] = competing_acceptor_key
					acceptor_val["competing_acceptor_dinucleotide"] = acceptor_dict[competing_acceptor_key]["dinucleotide"]
					acceptor_val["competing_acceptor_standard_dinucleotide"] = acceptor_dict[competing_acceptor_key]["standard_dinucleotide"]

					competitor_seq += [acceptor_dict[competing_acceptor_key]["seq"]]
					competitor_score += [acceptor_dict[competing_acceptor_key]["score"]]
					competitor_dinucleotide += [acceptor_dict[competing_acceptor_key]["dinucleotide"]]
					competitor_dinucleotide_standard += [acceptor_dict[competing_acceptor_key]["standard_dinucleotide"]]

				else:

					competitor_seq += ["NA"]
					competitor_score += ["NA"]
					competitor_dinucleotide += ["NA"]
					competitor_dinucleotide_standard += ["NA"]


				if donor_key in donor_dict:

					acceptor_val["donor_score"] = donor_dict[donor_key]["score"]
					acceptor_val["donor_seq"] = donor_dict[donor_key]["seq"]
					acceptor_val["donor_key"] = donor_key
					acceptor_val["donor_dinucleotide"] = donor_dict[donor_key]["dinucleotide"]
					acceptor_val["donor_standard_dinucleotide"] = donor_dict[donor_key]["standard_dinucleotide"]

					partner_seq += [donor_dict[donor_key]["seq"]]
					partner_score += [donor_dict[donor_key]["score"]]
					partner_dinucleotide += [donor_dict[donor_key]["dinucleotide"]]
					partner_dinucleotide_standard += [donor_dict[donor_key]["standard_dinucleotide"]]

				else:

					partner_seq += ["NA"]
					partner_score += ["NA"]
					partner_dinucleotide += ["NA"]
					partner_dinucleotide_standard += ["NA"]



		if "competitive_unique_donors" in event_val:

			for donor, donor_val in event_val["competitive_unique_donors"].iteritems():

				events += [event]

				competitive_junction_index += [donor_val["competitive_junction_index"]]

				donor_key = chrom + "_" + str(donor) + "_" + strand
				competing_donor_key = chrom + "_" + str(donor_val["competing_donor"]) + "_" + strand
				acceptor_key = chrom + "_" + str(donor_val["acceptor"]) + "_" + strand

				splice_site_type += ["donor"]
				primary_key += [donor_key]
				competitor_key += [competing_donor_key]
				partner_key += [acceptor_key]
				unique_form += [donor_val["unique_form"]]


				if donor_key in donor_dict:

					donor_val["donor_score"] = donor_dict[donor_key]["score"]
					donor_val["donor_seq"] = donor_dict[donor_key]["seq"]
					donor_val["donor_key"] = donor_key
					donor_val["donor_dinucleotide"] = donor_dict[donor_key]["dinucleotide"]
					donor_val["donor_standard_dinucleotide"] = donor_dict[donor_key]["standard_dinucleotide"]

					primary_seq += [donor_dict[donor_key]["seq"]]
					primary_score += [donor_dict[donor_key]["score"]]
					primary_dinucleotide += [donor_dict[donor_key]["dinucleotide"]]
					primary_dinucleotide_standard += [donor_dict[donor_key]["standard_dinucleotide"]]


				else:

					primary_seq += ["NA"]
					primary_score += ["NA"]
					primary_dinucleotide += ["NA"]
					primary_dinucleotide_standard += ["NA"]	


				if competing_donor_key in donor_dict:

					donor_val["competing_donor_score"] = donor_dict[competing_donor_key]["score"]
					donor_val["competing_donor_seq"] = donor_dict[competing_donor_key]["seq"]
					donor_val["competing_donor_key"] = competing_donor_key
					donor_val["competing_donor_dinucleotide"] = donor_dict[competing_donor_key]["dinucleotide"]
					donor_val["competing_donor_standard_dinucleotide"] = donor_dict[competing_donor_key]["standard_dinucleotide"]

					competitor_seq += [donor_dict[competing_donor_key]["seq"]]
					competitor_score += [donor_dict[competing_donor_key]["score"]]
					competitor_dinucleotide += [donor_dict[competing_donor_key]["dinucleotide"]]
					competitor_dinucleotide_standard += [donor_dict[competing_donor_key]["standard_dinucleotide"]]


				else:

					competitor_seq += ["NA"]
					competitor_score += ["NA"]
					competitor_dinucleotide += ["NA"]
					competitor_dinucleotide_standard += ["NA"]


				if acceptor_key in acceptor_dict:

					donor_val["acceptor_score"] = acceptor_dict[acceptor_key]["score"]
					donor_val["acceptor_seq"] = acceptor_dict[acceptor_key]["seq"]
					donor_val["acceptor_key"] = acceptor_key
					donor_val["acceptor_dinucleotide"] = acceptor_dict[acceptor_key]["dinucleotide"]
					donor_val["acceptor_standard_dinucleotide"] = acceptor_dict[acceptor_key]["standard_dinucleotide"]

					partner_seq += [acceptor_dict[acceptor_key]["seq"]]
					partner_score += [acceptor_dict[acceptor_key]["score"]]
					partner_dinucleotide += [acceptor_dict[acceptor_key]["dinucleotide"]]
					partner_dinucleotide_standard += [acceptor_dict[acceptor_key]["standard_dinucleotide"]]


				else:

					partner_seq += ["NA"]
					partner_score += ["NA"]
					partner_dinucleotide += ["NA"]
					partner_dinucleotide_standard += ["NA"]


	splice_sites_df = pd.DataFrame({"event_id": events,
				  "unique_form": unique_form, 
				  "splice_site_type": splice_site_type,
				  "competitive_junction_index": competitive_junction_index,
				   species + "_primary_key": primary_key, 
				   species + "_primary_seq": primary_seq,
				   species + "_primary_score": primary_score,
				   species + "_primary_dinucleotide": primary_dinucleotide,
				   species + "_primary_dinucleotide_standard": primary_dinucleotide_standard,
				   species + "_competitor_key": competitor_key,
				   species + "_competitor_seq": competitor_seq,
				   species + "_competitor_score": competitor_score,
				   species + "_competitor_dinucleotide": competitor_dinucleotide,
				   species + "_competitor_dinucleotide_standard": competitor_dinucleotide_standard,
				   species + "_partner_key": partner_key,
				   species + "_partner_seq": partner_seq,
				   species + "_partner_score": partner_score,
				   species + "_partner_dinucleotide": partner_dinucleotide,
				   species + "_partner_dinucleotide_standard": partner_dinucleotide_standard})[["event_id",
																								"unique_form",
																								"splice_site_type",
																								"competitive_junction_index",
																								 species + "_primary_key",
																								 species + "_primary_seq",
																								 species + "_primary_score",
																								 species + "_primary_dinucleotide",
																								 species + "_primary_dinucleotide_standard",
																								 species + "_competitor_key",
																								 species + "_competitor_seq",
																								 species + "_competitor_score",
																								 species + "_competitor_dinucleotide",
																								 species + "_competitor_dinucleotide_standard",
																								 species + "_partner_key",
																								 species + "_partner_seq",
																								 species + "_partner_score",
																								 species + "_partner_dinucleotide",
																								 species + "_partner_dinucleotide_standard"]]

	splice_sites_df[
		species + 
		"_primary_score"
		] = pd.to_numeric(splice_sites_df[species + "_primary_score"], 
		                  errors = "coerce")

	splice_sites_df[
		species + 
		"_competitor_score"
		] = pd.to_numeric(splice_sites_df[species + "_competitor_score"], 
		                                  errors = "coerce")

	splice_sites_df[
		species + 
		"_partner_score"
		] = pd.to_numeric(splice_sites_df[species + "_partner_score"], 
		                  errors = "coerce")

	splice_sites_df[
		species + 
		"_competitor_delta"
		] = abs(splice_sites_df[species + "_competitor_score"] - 
		        splice_sites_df[species + "_primary_score"])

	outfile_path = outdir + "/" + species + "_splice_sites.tsv"

	with open(outfile_path, 'w') as file:

		file.write("# " + full_command_line + "\n")

	splice_sites_df.to_csv(
		outfile_path, 
		sep = "\t", 
		index = False,
		mode = "a")

	return splice_sites_df


def main():

	full_command_line = " ".join(sys.argv)

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--outdir", 
		type = str, 
		help = "Path to output directory.", 
		required = True)

	parser.add_argument(
		"--species", 
		type = str, 
		help = ("Comma separated list of species e.g. " + 
			    "--species human,chimpanzee,orangutan"), 
		required = True)

	parser.add_argument(
		"--event_gtfs", 
		type = str, 
		help = ("Path to splice_lib event gtfs, " + 
			    "comma separated in same order as species."), 
		required = True)

	parser.add_argument(
		"--genome_fastas", 
		type = str, 
		help = ("Path to genome fasta files, " + 
			    "comma separated in same order as species."), 
		required = True)

	parser.add_argument(
		"--region_bed",
		action = "store_true",
		help = ("If set, regions surrounding competing splice " + 
				"sites will be written to bed file.  Useful for " + 
				"generating metaplots"))

	parser.add_argument(
		"--exonic_length",
		type = int,
		help = ("Number of exonic nucleotides to include " + 
			    "in region_bed output. Default = 30"),
		default = 30)

	parser.add_argument(
		"--intronic_length",
		type = int,
		help = ("Number of intronic nucleotides to include " + 
			    "in region_bed output. Default = 150"),
		default = 150)

	parser.add_argument(
		"--bedtools_path", 
		type = str, 
		help = "Path to bedtools exec. Default = 'bedtools'", 
		default = "bedtools")


	args = parser.parse_args()

	species_list = args.species.split(",")
	genome_fastas = args.genome_fastas.split(",")
	event_gtfs = args.event_gtfs.split(",")

	species_comparisons = []

	for i in range(0, len(species_list)):

		if i < len(species_list) - 1:

			for j in range(i + 1, len(species_list)):

				species_comparisons += [[species_list[i], species_list[j]]]


	species_dfs = []


	for i, species in enumerate(species_list):

		donor_coordinates = set()
		acceptor_coordinates = set()

		event_dict = splice_lib.generate_standard_event_dict(event_gtfs[i])
		splice_lib.complete_event_dict(event_dict)

		for event_val in event_dict.itervalues():

			get_competing_splice_sites(
				event_val, 
				acceptor_coordinates, 
				donor_coordinates)


		write_splice_site_bed(
			donor_coordinates, 
			"donor", 
			args.outdir, 
			species)

		write_splice_site_bed(
			acceptor_coordinates, 
			"acceptor", 
			args.outdir, 
			species)


		getfasta(
			args.outdir, 
			genome_fastas[i], 
			species, 
			args.bedtools_path)

		(donor_dict, 
		 acceptor_dict) = read_and_score_fasta(
						 	args.outdir, 
						 	species)

		species_dfs += [
		associate_seqs_scores_with_events(
			event_dict, 
			donor_dict, 
			acceptor_dict, 
			args.outdir, 
			species,
			full_command_line)
		]

		if args.region_bed:

			write_splice_site_region_bed(
				event_dict,
				args.outdir,
				species,
				exonic = args.exonic_length,
				intronic = args.intronic_length)

	if len(species_list) > 1:

		all_species_df = reduce(
			lambda df1, df2: df1.merge(
				df2, 
			    how = "inner", 
			    on = ["event_id", 
			   		  "unique_form",
			          "splice_site_type", 
			          "competitive_junction_index"]), 
			    species_dfs)

		for comparison in species_comparisons:

			key_base = comparison[0] + "_" + comparison[1]

			all_species_df[
				key_base + 
				"_same_primary_dinucleotide"
				] = (all_species_df[comparison[0] + "_primary_dinucleotide"] == 
					 all_species_df[comparison[1] + "_primary_dinucleotide"])

			all_species_df[
				key_base + 
				"_same_competitor_dinucleotide"
				] = (all_species_df[comparison[0] + "_competitor_dinucleotide"] == 
					 all_species_df[comparison[1] + "_competitor_dinucleotide"])

			all_species_df[
				key_base + 
				"_same_partner_dinucleotide"
				] = (all_species_df[comparison[0] + "_partner_dinucleotide"] == 
					 all_species_df[comparison[1] + "_partner_dinucleotide"])

			all_species_df[
				key_base + 
				"_competitor_delta"
				] = (all_species_df[comparison[0] + "_competitor_delta"] - 
					 all_species_df[comparison[1] + "_competitor_delta"])


		outfile_path = args.outdir + "/all_species_splice_sites.tsv"

		with open(outfile_path, "w") as file:

			file.write("# " + full_command_line + "\n") 
			## writes command line call to beginning of output file

		all_species_df.to_csv(
			outfile_path, 
			sep = "\t", 
			index = False,
			mode = "a")


if __name__ == "__main__":

	main()



