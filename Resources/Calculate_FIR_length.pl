#!/usr/bin/perl

use strict;
use warnings;

#--------------------------------------------------------------------------
# CalculateFIR_length.pl
#--------------------------------------------------------------------------
#
# Calculates flanking intergenic regions for each gene per contig/chromosome/scaffold :
#
#Usage: perl CalculateFIR_length.pl
#
# Written by Diane Saunders and Joe Win, Kamoun Lab, The Sainsbury Laboratory, Norwich, UK.
# Date: 21/03/2013
#
#--------------------------------------------------------------------------


#Input file format. gff or gtf?

print "\n\*\* Is the file in gff3 or gtf format? (gff3/gtf):";
my $file_format = <STDIN>;
chomp $file_format;

#Specify name of gtf/gff file

print "\n\*\* Please enter the name of the gff or gtf file: ";
my $file = <STDIN>;


#Specify type or feature [must be gene, exon or start and stop codons]
print "\n\*\* Which feature/type to process? (gene/mRNA/exon):";
my $field_type = <STDIN>;
chomp $field_type;


#Specify output file name

print "\n\*\* Enter name of output file:";
my $outfile = <STDIN>;
open(OUTFILE, ">$outfile") or die;

my $first_time = "YES";
my $first_time_overlaps = "YES";
my $previous_contig = "";
my %exon_hash = ();
my %gff_info = ();
my %info = ();
my %contigs_for_FIR = ();
my %sorted_contigs = ();
my $left_distance = "";
my $right_distance = "";


#Process files in gff format. 

if ($file_format eq "gff3") {
	if ($field_type eq "gene") {
		process_gff_gene ($file, \%gff_info);
	} elsif ($field_type eq "mRNA") {
		process_gff_mRNA ($file, \%gff_info);
	} elsif ($field_type eq "exon") {
		process_gff_exon ($file, \%gff_info);
	} else {
		print "\n\*\* Don't understand feature/type entered!!!\n";
	}
	
#Process files in gtf format.

} elsif ($file_format eq "gtf") {
	if ($field_type eq "gene") {
		process_gtf_gene ($file, \%gff_info);
	} elsif ($field_type eq "mRNA") {
		process_gtf_mRNA ($file, \%gff_info);
	} elsif ($field_type eq "exon") {
		process_gtf_exon ($file, \%gff_info);
	} else {
		print "\n\*\* Don't understand feature/type entered!!!\n";
	}

#Print error if file format is not "gff" or "gtf"

} else {
	print "\n\*\* Don't understand file format entered!!!\n";
}

#Print header 

print OUTFILE "\"geneid\",\"strand\",\"fiveprime\",\"threeprime\"\n";

#Sort gene order along contig for printing

foreach my $contig (sort keys %gff_info) {
	sort_this_contig ($contig, \%gff_info, \%sorted_contigs);
	%contigs_for_FIR = %sorted_contigs;
	sort_overlaps ($contig, \%gff_info, \%sorted_contigs, \%contigs_for_FIR);
	calculate_FIR ($contig, \%contigs_for_FIR);

}
my $this_gene = "";
my $first = "";
foreach my $contig (sort keys %sorted_contigs) {
	foreach my $this_gene (keys %{$sorted_contigs{$contig}}) {
		print OUTFILE "\"$this_gene\"\,\"$contigs_for_FIR{$contig}{$this_gene}{'strand'}\"\,$contigs_for_FIR{$contig}{$this_gene}{'left'}\,$contigs_for_FIR{$contig}{$this_gene}{'right'}\n";
		#print OUTFILE "$contig\t$this_gene\t$sorted_contigs{$contig}{$this_gene}{'previous_gene'}\t$sorted_contigs{$contig}{$this_gene}{'next_gene'}\t$sorted_contigs{$contig}{$this_gene}{'start'}\t$sorted_contigs{$contig}{$this_gene}{'end'}\t$sorted_contigs{$contig}{$this_gene}{'strand'}\n";
		#print OUTFILE "$contig\t$this_gene\t$contigs_for_FIR{$contig}{$this_gene}{'start'}\t$contigs_for_FIR{$contig}{$this_gene}{'end'}\t$contigs_for_FIR{$contig}{$this_gene}{'strand'}\t$contigs_for_FIR{$contig}{$this_gene}{'left'}\t$contigs_for_FIR{$contig}{$this_gene}{'right'}\n";
	}
}
close (OUTFILE);
exit;

#--------------------------------------------------------------------------
# Subroutines
#--------------------------------------------------------------------------

#Order genes along each contig

sub sort_this_contig {
	my ($contig_id, $gff_info_ref, $sorted_contigs_ref) = @_;
	foreach my $id (keys %{$$gff_info_ref{$contig_id}}) {
		if (exists $$sorted_contigs_ref{$contig_id}) {
			my $first_gene = "";
			foreach my $this_id (keys %{$$sorted_contigs_ref{$contig_id}}) {
				if ($$sorted_contigs_ref{$contig_id}{$this_id}{'first'} eq "YES") {
					$first_gene = $this_id;
					last;
				}
			}
			my $current_gene = $first_gene;
			my $found = "NO";
			while ($found eq "NO") {
				if ($$gff_info_ref{$contig_id}{$id}{'start'} < $$sorted_contigs_ref{$contig_id}{$current_gene}{'start'}) {
					if ($$sorted_contigs_ref{$contig_id}{$current_gene}{'first'} eq "YES") {
						$$sorted_contigs_ref{$contig_id}{$current_gene}{'first'} = "NO";
						$$sorted_contigs_ref{$contig_id}{$current_gene}{'previous_gene'} = $id;
						$$sorted_contigs_ref{$contig_id}{$id}{'last'} = "NO";
						$$sorted_contigs_ref{$contig_id}{$id}{'first'} = "YES";
						$$sorted_contigs_ref{$contig_id}{$id}{'start'} = $$gff_info_ref{$contig_id}{$id}{'start'};
						$$sorted_contigs_ref{$contig_id}{$id}{'end'} = $$gff_info_ref{$contig_id}{$id}{'end'};
						$$sorted_contigs_ref{$contig_id}{$id}{'strand'} = $$gff_info_ref{$contig_id}{$id}{'strand'};
						$$sorted_contigs_ref{$contig_id}{$id}{'previous_gene'} = "NIL";
						$$sorted_contigs_ref{$contig_id}{$id}{'next_gene'} = $current_gene;
						last;
					} else {
						$current_gene = $$sorted_contigs_ref{$contig_id}{$current_gene}{'previous_gene'}; 
						next;
					}
				} elsif ($$gff_info_ref{$contig_id}{$id}{'start'} >= $$sorted_contigs_ref{$contig_id}{$current_gene}{'start'}) {
					if ($$sorted_contigs_ref{$contig_id}{$current_gene}{'last'} eq "YES"){
						$$sorted_contigs_ref{$contig_id}{$current_gene}{'last'} = "NO";
						$$sorted_contigs_ref{$contig_id}{$current_gene}{'next_gene'} = $id;
						$$sorted_contigs_ref{$contig_id}{$id}{'first'} = "NO";
						$$sorted_contigs_ref{$contig_id}{$id}{'start'} = $$gff_info_ref{$contig_id}{$id}{'start'};
						$$sorted_contigs_ref{$contig_id}{$id}{'end'} = $$gff_info_ref{$contig_id}{$id}{'end'};
						$$sorted_contigs_ref{$contig_id}{$id}{'strand'} = $$gff_info_ref{$contig_id}{$id}{'strand'};
						$$sorted_contigs_ref{$contig_id}{$id}{'previous_gene'} = $current_gene;
						$$sorted_contigs_ref{$contig_id}{$id}{'last'} = "YES";
						$$sorted_contigs_ref{$contig_id}{$id}{'next_gene'} = "NIL";
						last;
					} else {
						$current_gene = $$sorted_contigs_ref{$contig_id}{$current_gene}{'next_gene'};
						if ($$gff_info_ref{$contig_id}{$id}{'start'} <= $$sorted_contigs_ref{$contig_id}{$current_gene}{'start'}) {
							$found = "YES";
							my $prev_gene = $$sorted_contigs_ref{$contig_id}{$current_gene}{'previous_gene'};
							$$sorted_contigs_ref{$contig_id}{$prev_gene}{'next_gene'} = $id;
							$$sorted_contigs_ref{$contig_id}{$current_gene}{'previous_gene'} = $id;
							$$sorted_contigs_ref{$contig_id}{$id}{'first'} = "NO";
							$$sorted_contigs_ref{$contig_id}{$id}{'last'} = "NO";
							$$sorted_contigs_ref{$contig_id}{$id}{'previous_gene'} = $prev_gene;
							$$sorted_contigs_ref{$contig_id}{$id}{'next_gene'} = $current_gene;
							$$sorted_contigs_ref{$contig_id}{$id}{'start'} = $$gff_info_ref{$contig_id}{$id}{'start'};
							$$sorted_contigs_ref{$contig_id}{$id}{'end'} = $$gff_info_ref{$contig_id}{$id}{'end'};
							$$sorted_contigs_ref{$contig_id}{$id}{'strand'} = $$gff_info_ref{$contig_id}{$id}{'strand'};
						}
					}
				}
			}
		} else {
			$$sorted_contigs_ref{$contig_id}{$id}{'first'} = "YES";
			$$sorted_contigs_ref{$contig_id}{$id}{'start'} = $$gff_info_ref{$contig_id}{$id}{'start'};
			$$sorted_contigs_ref{$contig_id}{$id}{'end'} = $$gff_info_ref{$contig_id}{$id}{'end'};
			$$sorted_contigs_ref{$contig_id}{$id}{'strand'} = $$gff_info_ref{$contig_id}{$id}{'strand'};
			$$sorted_contigs_ref{$contig_id}{$id}{'next_gene'} = "NIL";
			$$sorted_contigs_ref{$contig_id}{$id}{'previous_gene'} = "NIL";
			$$sorted_contigs_ref{$contig_id}{$id}{'last'} = "YES";
		}
	}
}
#--------------------------------------------------------------------------

#Calculate 5FIR and 3FIR distances

sub calculate_FIR {
	my ($contig_id, $contigs_for_FIR_ref) = @_;
	my $prev_id = "";
	my $next_id = "";
	my %overlap_genes = ();
	foreach my $id (keys %{$$contigs_for_FIR_ref{$contig_id}}) {
		if ($$contigs_for_FIR_ref{$contig_id}{$id}{'first'} eq "YES") {
			if ($$contigs_for_FIR_ref{$contig_id}{$id}{'next_gene'} eq "NIL") {
				$$contigs_for_FIR_ref{$contig_id}{$id}{'left'} = "NA";
				$$contigs_for_FIR_ref{$contig_id}{$id}{'right'} = "NA";
			} else {
				$next_id = $$contigs_for_FIR_ref{$contig_id}{$id}{'next_gene'};
				if ($$contigs_for_FIR_ref{$contig_id}{$id}{'strand'} eq "\+") {
					$$contigs_for_FIR_ref{$contig_id}{$id}{'left'} = "NA";
					$$contigs_for_FIR_ref{$contig_id}{$id}{'right'} = (($$contigs_for_FIR_ref{$contig_id}{$next_id}{'start'} - $$contigs_for_FIR_ref{$contig_id}{$id}{'end'}) - 1);
				} else {
					$$contigs_for_FIR_ref{$contig_id}{$id}{'right'} = "NA";
					$$contigs_for_FIR_ref{$contig_id}{$id}{'left'} = (($$contigs_for_FIR_ref{$contig_id}{$next_id}{'start'} - $$contigs_for_FIR_ref{$contig_id}{$id}{'end'}) -1);
				}	
			}
		} elsif ($$contigs_for_FIR_ref{$contig_id}{$id}{'last'} eq "YES") {
			if ($$contigs_for_FIR_ref{$contig_id}{$id}{'previous_gene'} eq "NIL") {
				$$contigs_for_FIR_ref{$contig_id}{$id}{'left'} = "NA";
				$$contigs_for_FIR_ref{$contig_id}{$id}{'right'} = "NA";
			} else {
				$prev_id = $$contigs_for_FIR_ref{$contig_id}{$id}{'previous_gene'};
				if ($$contigs_for_FIR_ref{$contig_id}{$id}{'strand'} eq "\+") {
					$$contigs_for_FIR_ref{$contig_id}{$id}{'left'} = (($$contigs_for_FIR_ref{$contig_id}{$id}{'start'} - $$contigs_for_FIR_ref{$contig_id}{$prev_id}{'end'}) -1);
					$$contigs_for_FIR_ref{$contig_id}{$id}{'right'} = "NA";	
				} else {
					$$contigs_for_FIR_ref{$contig_id}{$id}{'left'} = "NA";
					$$contigs_for_FIR_ref{$contig_id}{$id}{'right'} = (($$contigs_for_FIR_ref{$contig_id}{$id}{'start'} - $$contigs_for_FIR_ref{$contig_id}{$prev_id}{'end'}) -1);
				}
			}
		} elsif ($$contigs_for_FIR_ref{$contig_id}{$id}{'first'} eq "NO" && $$contigs_for_FIR_ref{$contig_id}{$id}{'last'} eq "NO") {
			$prev_id = $$contigs_for_FIR_ref{$contig_id}{$id}{'previous_gene'};
			$next_id = $$contigs_for_FIR_ref{$contig_id}{$id}{'next_gene'};
			if ($$contigs_for_FIR_ref{$contig_id}{$id}{'strand'} eq "\+") {
				$$contigs_for_FIR_ref{$contig_id}{$id}{'left'} = (($$contigs_for_FIR_ref{$contig_id}{$id}{'start'} - $$contigs_for_FIR_ref{$contig_id}{$prev_id}{'end'}) -1);
				$$contigs_for_FIR_ref{$contig_id}{$id}{'right'} = (($$contigs_for_FIR_ref{$contig_id}{$next_id}{'start'} - $$contigs_for_FIR_ref{$contig_id}{$id}{'end'}) -1);
			} else {
				$$contigs_for_FIR_ref{$contig_id}{$id}{'right'} = (($$contigs_for_FIR_ref{$contig_id}{$id}{'start'} - $$contigs_for_FIR_ref{$contig_id}{$prev_id}{'end'}) -1);
				$$contigs_for_FIR_ref{$contig_id}{$id}{'left'} = (($$contigs_for_FIR_ref{$contig_id}{$next_id}{'start'} - $$contigs_for_FIR_ref{$contig_id}{$id}{'end'}) -1);
			}

		}
	}
}

#--------------------------------------------------------------------------

#Deal with genes that overlap

sub sort_overlaps {
	my ($contig_id, $gff_info_ref, $sorted_contigs_ref, $contigs_for_FIR_ref) = @_;
	my $first_gene = "";
	
	# Finds the first gene on the contig
	foreach my $geness (keys %{$$sorted_contigs_ref{$contig_id}}) {
		if ($$sorted_contigs_ref{$contig_id}{$geness}{'first'} eq "YES") {
			$first_gene = $geness;
			if (($$sorted_contigs_ref{$contig_id}{$first_gene}{'next_gene'} eq "NIL") || ($$sorted_contigs_ref{$contig_id}{$first_gene}{'last'} eq "YES")) {
				return;
			}
			last;
		}
	}
	
	# Process rest of the genes on the contig
	my $curr_gene = $first_gene;
	my $done = 0;
	my $overlap = "";
	while (!$done) {
		my $gene_before = $$sorted_contigs_ref{$contig_id}{$curr_gene}{'previous_gene'};
		if ($gene_before eq "NIL") {
			$curr_gene = $$sorted_contigs_ref{$contig_id}{$curr_gene}{'next_gene'};
			next;
		}
		# Checks if the current gene's start is within the previous gene which indicate overlap in gene
		if (($$sorted_contigs_ref{$contig_id}{$curr_gene}{'start'} >= $$sorted_contigs_ref{$contig_id}{$gene_before}{'start'}) &&
		($$sorted_contigs_ref{$contig_id}{$curr_gene}{'start'} <= $$sorted_contigs_ref{$contig_id}{$gene_before}{'end'})) {
			if ($overlap eq "") {
				$overlap = $gene_before."\t".$curr_gene;
			} else {
				$overlap .= "\t".$curr_gene;
			}		
			# Copies previous gene's start to current gene's start
			$$contigs_for_FIR_ref{$contig_id}{$curr_gene}{'start'} = $$sorted_contigs_ref{$contig_id}{$gene_before}{'start'};
			$$contigs_for_FIR_ref{$contig_id}{$curr_gene}{'first'} = $$sorted_contigs_ref{$contig_id}{$gene_before}{'first'};
			# Copies current gene's "next_gene" to previous gene's "next_gene"
			$$contigs_for_FIR_ref{$contig_id}{$gene_before}{'next_gene'} = $$sorted_contigs_ref{$contig_id}{$curr_gene}{'next_gene'};
			$$contigs_for_FIR_ref{$contig_id}{$gene_before}{'last'} = $$sorted_contigs_ref{$contig_id}{$curr_gene}{'last'};
			$$contigs_for_FIR_ref{$contig_id}{$curr_gene}{'previous_gene'} = $$sorted_contigs_ref{$contig_id}{$gene_before}{'previous_gene'};
			
			if ($$sorted_contigs_ref{$contig_id}{$curr_gene}{'end'} > $$sorted_contigs_ref{$contig_id}{$gene_before}{'end'}) {
				# Current gene is longer than the previous gene that overlaps with
				$$contigs_for_FIR_ref{$contig_id}{$gene_before}{'end'} = $$sorted_contigs_ref{$contig_id}{$curr_gene}{'end'};
				$$contigs_for_FIR_ref{$contig_id}{$gene_before}{'next_gene'} = $$sorted_contigs_ref{$contig_id}{$curr_gene}{'next_gene'};
				$$contigs_for_FIR_ref{$contig_id}{$gene_before}{'last'} = $$sorted_contigs_ref{$contig_id}{$curr_gene}{'last'};
			} else {
				# Current gene is located within the previous gene
				$$contigs_for_FIR_ref{$contig_id}{$curr_gene}{'end'} = $$sorted_contigs_ref{$contig_id}{$gene_before}{'end'};
				$$contigs_for_FIR_ref{$contig_id}{$curr_gene}{'next_gene'} = $$sorted_contigs_ref{$contig_id}{$gene_before}{'next_gene'};
				$$contigs_for_FIR_ref{$contig_id}{$curr_gene}{'last'} = $$sorted_contigs_ref{$contig_id}{$gene_before}{'last'};
			}
			if ($$sorted_contigs_ref{$contig_id}{$curr_gene}{'last'} eq "YES") {
				if ($overlap ne "") {
					#print "$overlap\n";
					fill_overlap ($overlap, $contig_id, $sorted_contigs_ref, $contigs_for_FIR_ref);
					$overlap = "";
				}
				$done = 1;
			} else {
				$curr_gene = $$sorted_contigs_ref{$contig_id}{$curr_gene}{'next_gene'};
			}
		} else {
			if ($overlap ne "") {
				#print "$overlap\n";
				fill_overlap ($overlap, $contig_id, $sorted_contigs_ref, $contigs_for_FIR_ref);
				$overlap = "";
			}
			if ($$sorted_contigs_ref{$contig_id}{$curr_gene}{'last'} eq "YES") {
				$done = 1;
			} else {
				$curr_gene = $$sorted_contigs_ref{$contig_id}{$curr_gene}{'next_gene'};
			}
		}
	}
}

#--------------------------------------------------------------------------

sub fill_overlap {
	my ($overlap, $contig_id, $sorted_contigs_ref, $contigs_for_FIR_ref) = @_;
	my @gene_ids = split "\t", $overlap;
	my $last_gene = pop (@gene_ids);

	foreach my $gene(@gene_ids) {
		$$contigs_for_FIR_ref{$contig_id}{$gene}{'next_gene'} = $$contigs_for_FIR_ref{$contig_id}{$last_gene}{'next_gene'};
		$$contigs_for_FIR_ref{$contig_id}{$gene}{'previous_gene'} = $$contigs_for_FIR_ref{$contig_id}{$last_gene}{'previous_gene'};
		$$contigs_for_FIR_ref{$contig_id}{$gene}{'first'} = $$contigs_for_FIR_ref{$contig_id}{$last_gene}{'first'};
		$$contigs_for_FIR_ref{$contig_id}{$gene}{'last'} = $$contigs_for_FIR_ref{$contig_id}{$last_gene}{'last'};
		$$contigs_for_FIR_ref{$contig_id}{$gene}{'start'} = $$contigs_for_FIR_ref{$contig_id}{$last_gene}{'start'};
		$$contigs_for_FIR_ref{$contig_id}{$gene}{'end'} = $$contigs_for_FIR_ref{$contig_id}{$last_gene}{'end'};			
		
	}
}

#--------------------------------------------------------------------------

sub process_gff_gene {
	my ($in_file, $gff_info_ref) = @_;
	open(FILE, "<$in_file") or die;
	while (my $line = <FILE>) {
		chomp $line;
		next if ($line =~ /^\#.*$/);
		my ($contig, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split ("\t", $line);
		if ($type eq "gene") {
			my ($ID) = $attributes =~ /ID=(.*)\;?.*$/;
			if ($start < $end) {
				$$gff_info_ref{$contig}{$ID}{'start'} = $start;
				$$gff_info_ref{$contig}{$ID}{'end'} = $end;
				$$gff_info_ref{$contig}{$ID}{'strand'} = $strand;
			} else {
				$$gff_info_ref{$contig}{$ID}{'start'} = $end;
				$$gff_info_ref{$contig}{$ID}{'end'} = $start;
				$$gff_info_ref{$contig}{$ID}{'strand'} = $strand;	
			}		
		}
	}
	close FILE;
}

#--------------------------------------------------------------------------

sub process_gff_mRNA {
my ($in_file, $gff_info_ref) = @_;
	open(FILE, "<$in_file") or die;
	while (my $line = <FILE>) {
		chomp $line;
		next if ($line =~ /^\#.*$/);
		my ($contig, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split ("\t", $line);
		if ($type eq "mRNA") {
			my ($ID) = $attributes =~ /ID=(.*)\;?.*$/;
			if ($start < $end) {
				$$gff_info_ref{$contig}{$ID}{'start'} = $start;
				$$gff_info_ref{$contig}{$ID}{'end'} = $end;
				$$gff_info_ref{$contig}{$ID}{'strand'} = $strand;
			} else {
				$$gff_info_ref{$contig}{$ID}{'start'} = $end;
				$$gff_info_ref{$contig}{$ID}{'end'} = $start;
				$$gff_info_ref{$contig}{$ID}{'strand'} = $strand;	
			}		
		}
	}
	close FILE;
}


#--------------------------------------------------------------------------

sub process_gff_exon {
	my ($in_file, $gff_info_ref) = @_;
	my %exon_hash;
		open(FILE, "<$in_file") or die;
		while (my $line = <FILE>) {
			chomp $line;
			next if ($line =~ /^\#.*$/); 
			my ($contig, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split ("\t", $line);
			if ($type eq "exon") {
			my ($ID) = $attributes =~ /ID=(.*)\;?.*$/;
			$exon_hash{$ID}{'contig'} = $contig;
			if ($strand eq "+") {
				$exon_hash{$ID}{'strand'} = $strand;
				$exon_hash{$ID}{'start'}{$start} = 1;
				$exon_hash{$ID}{'end'}{$end} = 1;
			}elsif ($strand eq "-") {
				$exon_hash{$ID}{'strand'} = $strand;
				$exon_hash{$ID}{'start'}{$end} = 1;
				$exon_hash{$ID}{'end'}{$start} = 1;
			}
		}
		
	}
	foreach my $id (sort {$a cmp $b} keys %exon_hash) {
			#print "Processing ID: $id\n";
			if ($exon_hash{$id}{'strand'} eq "\+") {
				my $start_exon = (sort {$a <=> $b} keys %{$exon_hash{$id}{'start'}})[0];
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'} = $start_exon;
				my $end_exon = (sort {$b <=> $a} keys %{$exon_hash{$id}{'end'}})[0];
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'} = $end_exon;
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'strand'} = $exon_hash{$id}{'strand'};
			} else {
				my $start_exon = (sort {$a <=> $b} keys %{$exon_hash{$id}{'end'}})[0];
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'} = $start_exon;
				my $end_exon = (sort {$a <=> $b} keys %{$exon_hash{$id}{'start'}})[0];
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'} = $end_exon;
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'strand'} = $exon_hash{$id}{'strand'};
			}
			if 	($$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'} < $$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'}) {
				my $temp_end = $$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'};
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'} = $$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'};
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'} = $temp_end;
				
			}
	}
	close FILE;
}

#--------------------------------------------------------------------------

sub process_gtf_gene {
	my ($in_file, $gff_info_ref) = @_;
	open(FILE, "<$in_file") or die;
	while (my $line = <FILE>) {
		chomp $line;
		next if ($line =~ /^\#.*$/);
		my ($contig, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split ("\t", $line);
		if ($type eq "gene") {
			my ($ID) = $attributes =~ /\S+ \"(.*)\"\; .*$/;
			if ($start < $end) {
				$$gff_info_ref{$contig}{$ID}{'start'} = $start;
				$$gff_info_ref{$contig}{$ID}{'end'} = $end;
				$$gff_info_ref{$contig}{$ID}{'strand'} = $strand;
			} else {
				$$gff_info_ref{$contig}{$ID}{'start'} = $end;
				$$gff_info_ref{$contig}{$ID}{'end'} = $start;
				$$gff_info_ref{$contig}{$ID}{'strand'} = $strand;	
			}		

		}
	}
}

#--------------------------------------------------------------------------

sub process_gtf_mRNA {
my ($in_file, $gff_info_ref) = @_;
	open(FILE, "<$in_file") or die;
	while (my $line = <FILE>) {
		chomp $line;
		next if ($line =~ /^\#.*$/);
		my ($contig, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split ("\t", $line);
		if ($type eq "mRNA") {
			my ($ID) = $attributes =~ /\S+ \"(.*)\"\; .*$/;
			if ($start < $end) {
				$$gff_info_ref{$contig}{$ID}{'start'} = $start;
				$$gff_info_ref{$contig}{$ID}{'end'} = $end;
				$$gff_info_ref{$contig}{$ID}{'strand'} = $strand;
			} else {
				$$gff_info_ref{$contig}{$ID}{'start'} = $end;
				$$gff_info_ref{$contig}{$ID}{'end'} = $start;
				$$gff_info_ref{$contig}{$ID}{'strand'} = $strand;	
			}		
		}
	}
	close FILE;
}

#--------------------------------------------------------------------------

sub process_gtf_exon {
	my ($in_file, $gff_info_ref) = @_;
	my %exon_hash;
		open(FILE, "<$in_file") or die;
		while (my $line = <FILE>) {
			chomp $line;
			next if ($line =~ /^\#.*$/); 
			my ($contig, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split ("\t", $line);
			if ($type eq "exon") {
			my ($ID) = $attributes =~ /\S+ \"(.*)\"\; .*$/;
			$exon_hash{$ID}{'contig'} = $contig;
			if ($strand eq "+") {
				$exon_hash{$ID}{'strand'} = $strand;
				$exon_hash{$ID}{'start'}{$start} = 1;
				$exon_hash{$ID}{'end'}{$end} = 1;
			}elsif ($strand eq "-") {
				$exon_hash{$ID}{'strand'} = $strand;
				$exon_hash{$ID}{'start'}{$end} = 1;
				$exon_hash{$ID}{'end'}{$start} = 1;
			}
		}
		
	}
	foreach my $id (sort {$a cmp $b} keys %exon_hash) {
			if ($exon_hash{$id}{'strand'} eq "\+") {
				my $start_exon = (sort {$a <=> $b} keys %{$exon_hash{$id}{'start'}})[0];
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'} = $start_exon;
				my $end_exon = (sort {$b <=> $a} keys %{$exon_hash{$id}{'end'}})[0];
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'} = $end_exon;
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'strand'} = $exon_hash{$id}{'strand'};
			} else {
				my $start_exon = (sort {$a <=> $b} keys %{$exon_hash{$id}{'end'}})[0];
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'} = $start_exon;
				my $end_exon = (sort {$b <=> $a} keys %{$exon_hash{$id}{'start'}})[0];
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'} = $end_exon;
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'strand'} = $exon_hash{$id}{'strand'};
			}	
			if 	($$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'} < $$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'}) {
				my $temp_end = $$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'};
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'end'} = $$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'};
				$$gff_info_ref{"$exon_hash{$id}{'contig'}"}{$id}{'start'} = $temp_end;
		}
	}
	close FILE;
}