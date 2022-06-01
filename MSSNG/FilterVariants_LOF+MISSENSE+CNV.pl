#!/usr/bin/perl

# Version 1
# Author : Jaqueline Wang
# Bioinformatician at Human Genome and Stem Cell Research Center - IB USP
# e-mail : jaqueytw@gmail.com
# Date : 11 February 2021

### SCRIPT TO: filter inherited variants

use strict;

my ($file1, $file2, $file3) = @ARGV; 

if (!($file1 && $file2 && $file3)) {
        die "Missing input \
Usage : FilterVariants_LOF+MISSENSE+CNV.pl <MISSENSE_file> <LOF_file> <CNV_file> \n";
}

### THE CODE


## Missense FILE
open (IN, "$file1") or die "Fail to open $file1\n";

# Get columns names from missense file
my $head = <IN>;
chomp ($head);
my @HEAD = split (/\t/, $head);

my $fun = 0; ++$fun until $HEAD[$fun] =~ m/^Func\.refGene$/ or $fun > $#HEAD;
my $exo = 0; ++$exo until $HEAD[$exo] =~ m/^ExonicFunc\.refGene$/ or $exo > $#HEAD;
my $gn1 = 0; ++$gn1 until $HEAD[$gn1] =~ m/^gnomAD_exome_ALL$/ or $gn1 > $#HEAD;
my $gn2 = 0; ++$gn2 until $HEAD[$gn2] =~ m/^gnomAD_genome_ALL$/ or $gn2 > $#HEAD;
my $kg1 = 0; ++$kg1 until $HEAD[$kg1] =~ m/^1000g2015aug_all$/ or $kg1 > $#HEAD;
my $abr = 0; ++$abr until $HEAD[$abr] =~ m/^60mais609_AC$/ or $abr > $#HEAD;
my $cad = 0; ++$cad until $HEAD[$cad] =~ m/^CADD_phred_MSSNG$/ or $cad > $#HEAD;
my $sam = 0; ++$sam until $HEAD[$sam] =~ m/^Individual\.ID$/ or $sam > $#HEAD;
my $fam = 0; ++$fam until $HEAD[$fam] =~ m/^Family\.ID$/ or $fam > $#HEAD;
my $gen = 0; ++$gen until $HEAD[$gen] =~ m/^gene_symbol$/ or $gen > $#HEAD;
my $rel = 0; ++$rel until $HEAD[$rel] =~ m/^Relation$/ or $rel > $#HEAD;
my $mom = 0; ++$mom until $HEAD[$mom] =~ m/^Mother\.ID$/ or $mom > $#HEAD;
my $dad = 0; ++$dad until $HEAD[$dad] =~ m/^Father\.ID$/ or $dad > $#HEAD;
my $sid = 0; ++$sid until $HEAD[$sid] =~ m/^SID$/ or $sid > $#HEAD;
my $sex = 0; ++$sex until $HEAD[$sex] =~ m/^Sex$/ or $sex > $#HEAD;

my %hash;
my %related;
while (my $line = <IN>) {
	chomp ($line);
	my @data = split (/\t/, $line);	

	if ($data[$rel] =~ m/(child|proband|sibling)/) {
		if (!exists($related{$data[$sam]})) {
			$related{$data[$sam]} = $data[$fam];
		}
		if (($data[$mom] ne "-") && (!exists($related{$data[$mom]}))) {
			$related{$data[$mom]} = $data[$fam];
		}
		if (($data[$dad] ne "-") && (!exists($related{$data[$dad]}))) {
			$related{$data[$dad]} = $data[$fam];
		}
	}
	else {
		if (!exists($related{$data[$sam]})) {
			$related{$data[$sam]} = $data[$fam];
		}
	}

	if ($data[$fun] =~ m/(exonic|splicing)/) {
		if ($data[$exo]	=~ m/(frameshift deletion|frameshift insertion|nonframeshift deletion|nonframeshift insertion|nonsynonymous SNV|stopgain|stoploss|NA)/) {
			if (($data[$gn1] <= 0.01) && ($data[$gn2] <= 0.01) && ($data[$kg1] <= 0.01) && ($data[$abr] <= 12) && ($data[$cad] >= 20)) {
				if (!exists($hash{$data[$fam]})) {
					$hash{$data[$fam]}{$data[$sid]}{ MIS } = "$line";
				}
				else {
					if (!exists($hash{$data[$fam]}{$data[$sid]})) {
						$hash{$data[$fam]}{$data[$sid]}{ MIS } = "$line";
                                        }
					else {
						$hash{$data[$fam]}{$data[$sid]}{ MIS } = "$hash{$data[$fam]}{$data[$sid]}{ MIS }\n$line";
					}
				}
			}
		}
	}
}
close (IN);


## CNV FIle
open (IN3, "$file3") or die "Fail to open $file3\n";

# Get columns names from cnv file
my $head_cnv = <IN3>;
chomp ($head_cnv);
my @HEAD_cnv = split (/\t/, $head_cnv);

my $sam_cnv = 0; ++$sam_cnv until $HEAD_cnv[$sam_cnv] =~ m/^SID$/ or $sam_cnv > $#HEAD_cnv;
my $gen_cnv = 0; ++$gen_cnv until $HEAD_cnv[$gen_cnv] =~ m/^gene\_symbol$/ or $gen_cnv > $#HEAD_cnv;
my $rel_cnv = 0; ++$rel_cnv until $HEAD_cnv[$rel_cnv] =~ m/^Relation$/ or $rel_cnv > $#HEAD_cnv;
my $fam_cnv = 0; ++$fam_cnv until $HEAD_cnv[$fam_cnv] =~ m/^Family.ID$/ or $fam_cnv > $#HEAD_cnv;
my $mom_cnv = 0; ++$mom_cnv until $HEAD_cnv[$mom_cnv] =~ m/^Mother.ID$/ or $mom_cnv > $#HEAD_cnv;
my $dad_cnv = 0; ++$dad_cnv until $HEAD_cnv[$dad_cnv] =~ m/^Father.ID$/ or $dad_cnv > $#HEAD_cnv;
my $sex_cnv = 0; ++$sex_cnv until $HEAD_cnv[$sex_cnv] =~ m/^Sex$/ or $sex_cnv > $#HEAD_cnv;


my @CNV_fail;

while (my $line = <IN3>) {
        chomp ($line);
        my @data = split (/\t/, $line);

	if ($data[$rel_cnv] =~ m/(child|proband|sibling)/) {
                if (!exists($related{$data[$sam_cnv]})) {
                        $related{$data[$sam_cnv]} = $data[$fam_cnv];
                }
                if (($data[$mom_cnv] ne "-") && (!exists($related{$data[$mom_cnv]}))) {
                        $related{$data[$mom]} = $data[$fam_cnv];
                }
                if (($data[$dad_cnv] ne "-") && (!exists($related{$data[$dad_cnv]}))) {
                        $related{$data[$dad_cnv]} = $data[$fam_cnv];
                }
        }
        else {
                if (!exists($related{$data[$sam_cnv]})) {
                        $related{$data[$sam_cnv]} = $data[$fam_cnv];
                }
        }

        if (exists($hash{$data[$fam_cnv]})) {
		if (!exists($hash{$data[$fam_cnv]}{$data[$sam_cnv]})) {
			$hash{$data[$fam_cnv]}{$data[$sam_cnv]}{ CNV } = "$line";
		}
		else {
			if (!exists($hash{$data[$fam_cnv]}{$data[$sam_cnv]}{ CNV })) {
				$hash{$data[$fam_cnv]}{$data[$sam_cnv]}{ CNV } = "$line";
			}
			else {
				$hash{$data[$fam_cnv]}{$data[$sam_cnv]}{ CNV } = "$hash{$data[$fam_cnv]}{$data[$sam_cnv]}{ CNV }\n$line";
			}
		}
        }
        else {
		$hash{$data[$fam_cnv]}{$data[$sam_cnv]}{ CNV } = "$line";
        }
}
close (IN3);


## LOF File
open (IN2, "$file2") or die "Fail to open $file2\n";

# Get columns names from missense file
my $head_lof = <IN2>;
chomp ($head_lof);
my @HEAD_lof = split (/\t/, $head_lof);

my $fun_lof = 0; ++$fun_lof until $HEAD_lof[$fun_lof] =~ m/^Func\.refGene$/ or $fun_lof > $#HEAD_lof;
my $exo_lof = 0; ++$exo_lof until $HEAD_lof[$exo_lof] =~ m/^ExonicFunc\.refGene$/ or $exo_lof > $#HEAD_lof;
my $gn1_lof = 0; ++$gn1_lof until $HEAD_lof[$gn1_lof] =~ m/^gnomAD_exome_ALL$/ or $gn1_lof > $#HEAD_lof;
my $gn2_lof = 0; ++$gn2_lof until $HEAD_lof[$gn2_lof] =~ m/^gnomAD_genome_ALL$/ or $gn2_lof > $#HEAD_lof;
my $kg1_lof = 0; ++$kg1_lof until $HEAD_lof[$kg1_lof] =~ m/^1000g2015aug_all$/ or $kg1_lof > $#HEAD_lof;
my $abr_lof = 0; ++$abr_lof until $HEAD_lof[$abr_lof] =~ m/^60mais609_AC$/ or $abr_lof > $#HEAD_lof;
my $sam_lof = 0; ++$sam_lof until $HEAD_lof[$sam_lof] =~ m/^sample$/ or $sam_lof > $#HEAD_lof;
my $gen_lof = 0; ++$gen_lof until $HEAD_lof[$gen_lof] =~ m/^gene_symbol$/ or $gen_lof > $#HEAD_lof;
my $rel_lof = 0; ++$rel_lof until $HEAD_lof[$rel_lof] =~ m/^Relation$/ or $rel_lof > $#HEAD_lof;
my $mom_lof = 0; ++$mom_lof until $HEAD_lof[$mom_lof] =~ m/^Mother.ID$/ or $mom_lof > $#HEAD_lof;
my $dad_lof = 0; ++$dad_lof until $HEAD_lof[$dad_lof] =~ m/^Father.ID$/ or $dad_lof > $#HEAD_lof;
my $sex_lof = 0; ++$sex_lof until $HEAD_lof[$sex_lof] =~ m/^Sex$/ or $sex_lof > $#HEAD_lof;

my @LOF_fail;

while (my $line = <IN2>) {
        chomp ($line);
        my @data = split (/\t/, $line);

        if ($data[$fun_lof] =~ m/(exonic|splicing)/) {
                if ($data[$exo_lof] =~ m/(frameshift deletion|frameshift insertion|nonframeshift deletion|nonframeshift insertion|nonsynonymous SNV|stopgain|stoploss|NA)/) {
                        if (($data[$gn1_lof] <= 0.01) && ($data[$gn2_lof] <= 0.01) && ($data[$kg1_lof] <= 0.01) && ($data[$abr_lof] <= 12)) {
			
				if (exists($related{$data[$sam_lof]})) {
					my $fam = $related{$data[$sam_lof]};
	
					if (exists($hash{$fam}{$data[$sam_lof]})) {
						if (!exists($hash{$fam}{$data[$sam_lof]}{ LOF })) {
							$hash{$fam}{$data[$sam_lof]}{ LOF } = "$line";
						}
						else {
							$hash{$fam}{$data[$sam_lof]}{ LOF } = "$hash{$fam}{$data[$sam_lof]}{ LOF }\n$line";
						}
					}
					else {
						$hash{$fam}{$data[$sam_lof]}{ LOF } = "$line";
					}
				}
				else {
					if ((exists($related{$data[$mom_lof]})) && ($mom_lof ne "-")) {
						my $fam = $related{$data[$mom_lof]};

						if (!exists($hash{$fam}{$data[$sam_lof]})) {
							$hash{$fam}{$data[$sam_lof]}{ LOF } = "$line";
						}
						else {
							if (!exists($hash{$fam}{$data[$sam_lof]}{ LOF })) {
								$hash{$fam}{$data[$sam_lof]}{ LOF } = "$line";
							}
							else {
								$hash{$fam}{$data[$sam_lof]}{ LOF } = "$hash{$fam}{$data[$sam_lof]}{ LOF }\n$line";
							}
						}
					}
					elsif ((exists($related{$data[$dad_lof]})) && ($dad_lof ne "-")) {
						my $fam = $related{$data[$dad_lof]};

						if (!exists($hash{$fam}{$data[$sam_lof]})) {
							$hash{$fam}{$data[$sam_lof]}{ LOF } = "$line";
						}
						else {
							if (!exists($hash{$fam}{$data[$sam_lof]}{ LOF })) {
								$hash{$fam}{$data[$sam_lof]}{ LOF } = "$line";
							}
							else {
								$hash{$fam}{$data[$sam_lof]}{ LOF } = "$hash{$fam}{$data[$sam_lof]}{ LOF }\n$line";
							}
						}
					}
					else {
						push @LOF_fail, $line;
					}


				}
			}
                }
        }
}
close (IN2);

# Genes via do Calcio
my %CAL;
open (IN, "Calcio.txt") or die "Fail to open calcio list\n";
while (my $line = <IN>) {
	chomp ($line);
	
	if (!exists($CAL{$line})) {
		$CAL{$line} = $line;
	}
	else {
		print "Repeated $line\n";
	}
}
close (IN);


# Genes via da Relina
my %REL;
open (IN, "Relina.txt") or die "Fail to open relina list\n";
while (my $line = <IN>) {
        chomp ($line);
        
        if (!exists($REL{$line})) {
                $REL{$line} = $line;
        }
        else {
                print "Repeated $line\n";
        }
}
close (IN);


foreach my $key1 (sort keys %hash) {

	my $cont = 0;
	foreach my $key2 (sort keys %{ $hash{$key1} }) {

		my $cal = 0; my $rel = 0;
		# Check Missense Variants
		if (exists($hash{$key1}{$key2}{ MIS })) {
			my @VAR = split (/\n/, $hash{$key1}{$key2}{ MIS });

			foreach my $var (@VAR) {
				my @data = split (/\t/, $var);

				if (exists($CAL{$data[$gen]})) {
					$cal += 1;
				}
				if (exists($REL{$data[$gen]})) {
					$rel += 1;
				}
				#print "$var\n";
			}
		}
		# Check Lof Variants
		if (exists($hash{$key1}{$key2}{ LOF })) {
			my @VAR = split (/\n/, $hash{$key1}{$key2}{ LOF });
			foreach my $var (@VAR) {
				my @data = split (/\t/, $var);

				if (exists($CAL{$data[$gen_lof]})) {
					$cal += 1;
				}
				if (exists($REL{$data[$gen_lof]})) {
					$rel += 1;
				}
			}
		}
		# Check CNV Variants
		if (exists($hash{$key1}{$key2}{ CNV })) {
			my @VAR = split (/\n/, $hash{$key1}{$key2}{ CNV });
			foreach my $var (@VAR) {
				my @data = split (/\t/, $var);

				if (exists($CAL{$data[$gen_cnv]})) {
					$cal += 1;
				}
				if (exists($REL{$data[$gen_cnv]})) {
					$rel += 1;
				}
			}
		}


		# If the Sample has at least 2 variants in REL list and 1 varians in CAL list
		if (($rel >= 2) && ($cal >= 1)) {
			$cont += 1;
		}
	}

	if ($cont >= 1) {
		print "FAMILY\tSAMPLE\tTYPE\tChr\tStart\tEnd\tRef\tAlt\tRelation\tGeneSymbol\tSex\n";

		foreach my $key2 (sort keys %{ $hash{$key1} }) {
			if (exists($hash{$key1}{$key2}{ MIS })) {
				#print "FAMILY\tSAMPLE\tTYPE\t$head\n";
				my @VAR = split (/\n/, $hash{$key1}{$key2}{ MIS });
				foreach my $var (@VAR) {
					my @data = split (/\t/, $var);
					print "$key1\t$key2\tMIS\t$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[$rel]\t$data[$gen]\t$data[$sex]\n";
				}
			}
			if (exists($hash{$key1}{$key2}{ LOF })) {
				#print "FAMILY\tSAMPLE\tTYPE\t$head_lof\n";
				my @VAR = split (/\n/, $hash{$key1}{$key2}{ LOF });
				foreach my $var (@VAR) {
					my @data = split (/\t/, $var);
					print "$key1\t$key2\tLOF\t$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[$rel_lof]\t$data[$gen_lof]\t$data[$sex_lof]\n";
				}
			}
			if (exists($hash{$key1}{$key2}{ CNV })) {
				#print "FAMILY\tSAMPLE\tTYPE\t$head_cnv\n";
				my @VAR = split (/\n/, $hash{$key1}{$key2}{ CNV });
				foreach my $var (@VAR) {
					my @data = split (/\t/, $var);
					print "$key1\t$key2\tCNV\t$data[2]\t$data[3]\t$data[4]\t$data[1]\t$data[5]\t$data[$rel_cnv]\t$data[$gen_cnv]\t$data[$sex_cnv]\n";
				}
			}
		}
		print "\n\n";
	}
}
