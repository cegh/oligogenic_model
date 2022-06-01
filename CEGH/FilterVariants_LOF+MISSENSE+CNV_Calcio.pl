#!/usr/bin/perl

# Version 1
# Author : Jaqueline Wang
# Bioinformatician at Human Genome and Stem Cell Research Center - IB USP
# e-mail : jaqueytw@gmail.com
# Date : 11 February 2021

### SCRIPT TO: filter inherited variants

use strict;

my ($vcf, $cnv1, $cnv2) = @ARGV; 

if (!($vcf && $cnv1 && $cnv2)) {
        die "Missing input \
Usage : FilterVariants_LOF+MISSENSE+CNV.pl <VCF_file> <CNV_females_file> <CNV_males_file> \n";
}

### THE CODE

# Realtionship file
my %relation;
open (IN, "/home/venus/mar/alunos/jaqueytw/ProfRita/Gabi/CEGH/Samples-Relationship_VCF.txt") or die "Fail to open relationship list\n";
while (my $line = <IN>) {
        chomp ($line);
	my @data = split (/\t/, $line);

        if (!exists($relation{$data[0]})) {
                $relation{$data[0]}{ fam } = $data[0];
		$relation{$data[0]}{ col } = 0;
        }
        else {
                print "Repeated $line\n";
        }

	if (!exists($relation{$data[1]})) {
		$relation{$data[1]}{ fam } = $data[0];
		$relation{$data[1]}{ col } = 0;
	}
	else {
		print "Repeated $line\n";
	}

	if (!exists($relation{$data[2]})) {
		$relation{$data[2]}{ fam } = $data[0];
		$relation{$data[2]}{ col } = 0;
	}
	else {
		print "Repeated $line\n";
	}

	if ($data[3] ne "OK") {
		if (!exists($relation{$data[3]})) {
			$relation{$data[3]}{ fam } = $data[0];
			$relation{$data[3]}{ col } = 0;
		}
		else {
			print "Repeated $line\n";
		}
	}
}
close (IN);


## VCF FILE
open (IN, "$vcf") or die "Fail to open $vcf\n";

# Get columns names from VCF file
my $head = <IN>;
chomp ($head);
my @HEAD = split (/\t/, $head);

my $fil = 0; ++$fil until $HEAD[$fil] =~ m/^FILTER$/ or $fil > $#HEAD;
my $for = 0; ++$for until $HEAD[$for] =~ m/^FORMAT$/ or $for > $#HEAD;
my $ge1 = 0; ++$ge1 until $HEAD[$ge1] =~ m/^Gene\.refGene$/ or $ge1 > $#HEAD;
my $ge2 = 0; ++$ge2 until $HEAD[$ge2] =~ m/^HGNC\_principal$/ or $ge2 > $#HEAD;
my $ge3 = 0; ++$ge3 until $HEAD[$ge3] =~ m/^HGNC\_previous$/ or $ge3 > $#HEAD;
my $fun = 0; ++$fun until $HEAD[$fun] =~ m/^Func\.refGene$/ or $fun > $#HEAD;
my $exo = 0; ++$exo until $HEAD[$exo] =~ m/^ExonicFunc\.refGene$/ or $exo > $#HEAD;
my $gn1 = 0; ++$gn1 until $HEAD[$gn1] =~ m/^gnomAD\_exome\_ALL$/ or $gn1 > $#HEAD;
my $gn2 = 0; ++$gn2 until $HEAD[$gn2] =~ m/^gnomAD\_genome\_ALL$/ or $gn2 > $#HEAD;
my $kg1 = 0; ++$kg1 until $HEAD[$kg1] =~ m/^1000g2015aug\_all$/ or $kg1 > $#HEAD;
my $abr = 0; ++$abr until $HEAD[$abr] =~ m/^60mais609\_AC$/ or $abr > $#HEAD;
my $cad = 0; ++$cad until $HEAD[$cad] =~ m/^CADD\_phred$/ or $cad > $#HEAD;
my $lof = 0; ++$lof until $HEAD[$lof] =~ m/^Loftee\_LoF$/ or $lof > $#HEAD;

my $col = $for + 1;
my %columns;
while ($col <= $#HEAD) {
	if (exists($relation{$HEAD[$col]})) {
		$relation{$HEAD[$col]}{ col } = $col;
	}
	else {
		print "Not exists $HEAD[$col]\n";
	}

	if (!exists($columns{$col})) {
		$columns{$col} = $HEAD[$col];
	}
	else {
		print "Repeated column $col\n";
	}

	$col += 1;
}


my %hash;
while (my $line = <IN>) {
	chomp ($line);
	my @data = split (/\t/, $line);	
	next unless $data[$fil] =~ m/^PASS$/;

	if ($data[$fun] =~ m/(exonic|splicing)/) {
		if (($data[$exo] =~ m/(frameshift deletion|frameshift insertion|nonframeshift deletion|nonframeshift insertion|nonsynonymous SNV|stopgain|stoploss|NA)/) || ($data[$lof] =~ m/(LC|HC){1}/)) {
			if (($data[$gn1] <= 0.01) && ($data[$gn2] <= 0.01) && ($data[$kg1] <= 0.01) && ($data[$abr] <= 12) && ($data[$cad] >= 20)) {
				my $cont = $for + 1;
				
				while ($cont <= $#data) {
					if ($data[$cont] =~ /^[01\.]{1}\/1:/) {
						my @new = @data[0 .. $for];
						push @new, $columns{$cont}, $data[$cont];
						my $ins = join ("\t", @new);

						#print "$relation{$columns{$cont}}{ fam }\t$columns{$cont}\t$data[0]\t$data[1]\t$data[2]\t$data[3]\n";

						if (!exists($hash{$relation{$columns{$cont}}{ fam }})) {
							$hash{$relation{$columns{$cont}}{ fam }}{$columns{$cont}}{ VCF } = "$ins";
						}
						else {
							if (!exists($hash{$relation{$columns{$cont}}{ fam }}{$columns{$cont}})) {
								$hash{$relation{$columns{$cont}}{ fam }}{$columns{$cont}}{ VCF } = "$ins";
							}
							else {
								$hash{$relation{$columns{$cont}}{ fam }}{$columns{$cont}}{ VCF } = "$hash{$relation{$columns{$cont}}{ fam }}{$columns{$cont}}{ VCF }\n$ins";
							}
						}
					}
					$cont += 1;
				}
			}
		}
	}
}
close (IN);


## CNV females
open (CNV1, "$cnv1") or die "Fail to open $cnv1\n";

# Get columns names from VCF file
my $headF = <CNV1>;
chomp ($headF);
my @CNV1 = split (/\t/, $headF);

my $cnv1_sam = 0; ++$cnv1_sam until $CNV1[$cnv1_sam] =~ m/^Samples\_ID$/ or $cnv1_sam > $#CNV1;
my $cnv1_for = 0; ++$cnv1_for until $CNV1[$cnv1_for] =~ m/^FORMAT$/ or $cnv1_for > $#CNV1;
my $cnv1_ann = 0; ++$cnv1_ann until $CNV1[$cnv1_ann] =~ m/^Annotation\_mode$/ or $cnv1_ann > $#CNV1;
my $cnv1_gen = 0; ++$cnv1_gen until $CNV1[$cnv1_gen] =~ m/^Gene_name$/ or $cnv1_gen > $#CNV1;

my %cnv1_cols;
my $col = $cnv1_for + 1;
while ($col < $cnv1_ann) {
	if (exists($relation{$CNV1[$col]})) {
                $relation{$CNV1[$col]}{ sex } = "F";
        }
        else {
                #print "Not exists $CNV1[$col]\n";
        }


	if (!exists($cnv1_cols{$CNV1[$col]})) {
		$cnv1_cols{$CNV1[$col]} = $col;
	}
	else {
		print "Repeated sample $CNV1[$col]\n";
	}
	$col += 1;
}

while (my $line = <CNV1>) {
	chomp ($line);
	my @data = split (/\t/, $line);
	next unless $data[$cnv1_ann] =~ m/^full$/;

	my @cnv_sam = split (/\,/, $data[$cnv1_sam]);

	if (scalar @cnv_sam <= 2) {
		foreach my $key (sort @cnv_sam) {
			my @new = @data[0 .. $cnv1_for];
			push @new, $data[$cnv1_cols{$key}], @data[$cnv1_ann .. $#data-1];
			my $ins = join ("\t", @new);

			if (exists($hash{$relation{$key}{ fam }}{$key})) {
				if (!exists($hash{$relation{$key}{ fam }}{$key}{ CNV })) {
					$hash{$relation{$key}{ fam }}{$key}{ CNV } = $ins;	
					#print "$ins\n";
				}
				else {
					$hash{$relation{$key}{ fam }}{$key}{ CNV } = "$hash{$relation{$key}{ fam }}{$key}{ CNV }\n$ins";
					#print "$ins\n";
				}
			}
			else {
				if (exists($hash{$relation{$key}{ fam }})) {
					if (!exists($hash{$relation{$key}{ fam }}{$key})) {
						$hash{$relation{$key}{ fam }}{$key}{ CNV } = $ins;
						#print "$key\t$data[$cnv1_for]\t$data[$cnv1_cols{$key}]\n";
					}		
				}
				else {
					print "$key\t$data[$cnv1_for]\t$data[$cnv1_cols{$key}]\n";
				}
			}	
		}
	}
}
close (CNV1);


## CNV males
open (CNV2, "$cnv2") or die "Fail to open $cnv2\n";

# Get columns names from VCF file
my $headM = <CNV2>;
chomp ($headM);
my @CNV2 = split (/\t/, $headM);

my $cnv2_sam = 0; ++$cnv2_sam until $CNV2[$cnv2_sam] =~ m/^Samples\_ID$/ or $cnv2_sam > $#CNV2;
my $cnv2_for = 0; ++$cnv2_for until $CNV2[$cnv2_for] =~ m/^FORMAT$/ or $cnv2_for > $#CNV2;
my $cnv2_ann = 0; ++$cnv2_ann until $CNV2[$cnv2_ann] =~ m/^Annotation\_mode$/ or $cnv2_ann > $#CNV2;
my $cnv2_gen = 0; ++$cnv2_gen until $CNV2[$cnv2_gen] =~ m/^Gene_name$/ or $cnv2_gen > $#CNV2;

my %cnv2_cols;
my $col = $cnv2_for + 1;
while ($col < $cnv2_ann) {
	if (exists($relation{$CNV2[$col]})) {
                $relation{$CNV2[$col]}{ sex } = "M";
        }
        else {
                #print "Not exists $CNV2[$col]\n";
        }

        if (!exists($cnv2_cols{$CNV2[$col]})) {
                $cnv2_cols{$CNV2[$col]} = $col;
        }
        else {
                print "Repeated sample $CNV2[$col]\n";
        }
        $col += 1;
}

while (my $line = <CNV2>) {
        chomp ($line);
        my @data = split (/\t/, $line);
        next unless $data[$cnv2_ann] =~ m/^full$/;

        my @cnv_sam = split (/\,/, $data[$cnv2_sam]);

        if (scalar @cnv_sam <= 2) {
                foreach my $key (sort @cnv_sam) {
                        my @new = @data[0 .. $cnv2_for];
                        push @new, $data[$cnv2_cols{$key}], @data[$cnv2_ann .. $#data-1];
                        my $ins = join ("\t", @new);

                        if (exists($hash{$relation{$key}{ fam }}{$key})) {
                                if (!exists($hash{$relation{$key}{ fam }}{$key}{ CNV })) {
                                        $hash{$relation{$key}{ fam }}{$key}{ CNV } = $ins;
                                        #print "$ins\n";
                                }
                                else {
                                        $hash{$relation{$key}{ fam }}{$key}{ CNV } = "$hash{$relation{$key}{ fam }}{$key}{ CNV }\n$ins";
                                        #print "$ins\n";
                                }
                        }
                        else {
                                if (exists($hash{$relation{$key}{ fam }})) {
                                        if (!exists($hash{$relation{$key}{ fam }}{$key})) {
                                                $hash{$relation{$key}{ fam }}{$key}{ CNV } = $ins;
                                                #print "$key\t$data[$cnv2_for]\t$data[$cnv2_cols{$key}]\n";
                                        }
                                }
                                else {
                                        print "$key\t$data[$cnv2_for]\t$data[$cnv2_cols{$key}]\n";
                                }
                        }
                }
        }
}
close (CNV2);


# Genes via do Calcio
my %CAL;
open (IN, "/home/venus/mar/alunos/jaqueytw/ProfRita/Gabi/AnnovarAtualizado/Concatenar/Calcio.txt") or die "Fail to open calcio list\n";
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
open (IN, "/home/venus/mar/alunos/jaqueytw/ProfRita/Gabi/AnnovarAtualizado/Concatenar/Relina.txt") or die "Fail to open relina list\n";
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


# Get CNV1 new columns
my @newCNV1 = @CNV1[0 .. $cnv1_for];
push @newCNV1, "Sample_info", @CNV1[$cnv1_ann .. $#CNV1-1];
my $newcnv1_gen = 0; ++$newcnv1_gen until $newCNV1[$newcnv1_gen] =~ m/^Gene_name$/ or $newcnv1_gen > $#newCNV1;
my $newcnv1_for = 0; ++$newcnv1_for until $newCNV1[$newcnv1_for] =~ m/^FORMAT$/ or $newcnv1_for > $#newCNV1;

# Get families that qualifies for the condition
foreach my $key1 (sort keys %hash) {

	my $cont = 0;
        foreach my $key2 (sort keys %{ $hash{$key1} }) {
		my $cal = 0; my $rel = 0;

		# Check VCF (Missense and LOF variants)
		my @VAR = split (/\n/, $hash{$key1}{$key2}{ VCF });
                foreach my $var (@VAR) {
                        my @data = split (/\t/, $var);
			my $cal1 = 0; my $rel1 = 0;
	
			my @GE1 = split (/;/, $data[$ge1]);
			foreach my $gen (@GE1) {
				if (exists($CAL{$gen})) {
					$cal1 += 1;
				}
				if (exists($REL{$gen})) {
					$rel1 += 1;
				}
			}

			my @GE2 = split (/\|/, $data[$ge2]);
			foreach my $gen (@GE2) {
				if (exists($CAL{$gen})) {
					$cal1 += 1;
				}
				if (exists($REL{$gen})) {
					$rel1 += 1;
				}
			}

			my @GE3 = split (/\|/, $data[$ge3]);
			foreach my $gen (@GE2) {
				if (exists($CAL{$gen})) {
					$cal1 += 1;
				}
				if (exists($REL{$gen})) {
					$rel1 += 1;
				}
			}

			if ($cal1 >= 1) {
				$cal += 1;
			}
			if ($rel1 >= 1) {
				$rel += 1;
			}
                }

		# Check CNV
		my @VAR = split (/\n/, $hash{$key1}{$key2}{ CNV });
                foreach my $var (@VAR) {
                        my @data = split (/\t/, $var);
                        my $cal1 = 0; my $rel1 = 0;

                        my @GE1 = split (/;/, $data[$newcnv1_gen]);
                        foreach my $gen (@GE1) {
                                if (exists($CAL{$gen})) {
                                        $cal1 += 1;
                                }
                                if (exists($REL{$gen})) {
                                        $rel1 += 1;
                                }
                        }

                        if ($cal1 >= 1) {
                                $cal += 1;
                        }
                        if ($rel1 >= 1) {
                                $rel += 1;
                        }
                }


		## If the Sample has at least 2 variants in REL list and 1 varians in CAL list
		#if (($rel >= 2) && ($cal >= 1)) {
		#	$cont += 1;
		#}

		# If the Sample has at least1 varians in CAL list
		if ($cal >= 1) {
			$cont += 1;
		}
        }

	if ($cont >= 1) {
		print "FAMILY\tSAMPLE\tTYPE\tChr\tStart\tEnd\tRef\tAlt\t\tFORMAT\tSample_INFO\tGeneSymbol\tSex\n";
		foreach my $key2 (sort keys %{ $hash{$key1} }) {
			my @VAR = split (/\n/, $hash{$key1}{$key2}{ VCF });
			foreach my $var (@VAR) {
				my @data = split (/\t/, $var);
				print "$key1\t$key2\tVCF\t$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[$for]\t$data[$for + 2]\t$data[$ge1]\t$relation{$key2}{ sex }\n";
			}

			my @VAR = split (/\n/, $hash{$key1}{$key2}{ CNV });
			foreach my $var (@VAR) {
                                my @data = split (/\t/, $var);
                                print "$key1\t$key2\tCNV\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[$newcnv1_for]\t$data[$newcnv1_for + 1]\t$data[$newcnv1_gen]\t$relation{$key2}{ sex }\n";
                        }
		}
		print "\n\n";
	}
}

