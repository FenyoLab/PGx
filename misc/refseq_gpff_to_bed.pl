use strict;
my $filename="";
my $filename2="";

my $error=0;

# [MANOR] Hardcoding done here:
$filename="proteome.gpff";
$filename2="refSeqAli.txt"; 

if ($error==0)
{
	my $line="";
	my %alignment_strand=();
	my %alignment_chr=();
	my %alignment_pos=();
	my %alignment_seg_len=();
	my %alignment_seg_pos=();
	my %description=();
	my %genes=();
	my %gis=();
	my %protein_transcript=();
	
	open (OUT,">$filename2.bed");
	if (open (IN,"$filename2"))
	{
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
			{
				$alignment_strand{$11}=$10;
				$alignment_chr{$11}=$15;
				$alignment_pos{$11}=$17;
				$alignment_seg_len{$11}="$20,";
				$alignment_seg_pos{$11}="$22,";
				my $start=$17;
				my $out="$15\t$17\t$18\t$11\t1000\t$10\t$17\t$18\t0\t$19";
				my $pos="";
				my $pos_ori="$22,";
				my $len="";
				my $len_ori="$20,";
				while($pos_ori=~s/^([^\,]+)\,//)
				{
					my $pos_=$1-$start;
					if ($pos=~/\w/) { $pos.=","; }
					$pos.="$pos_";
				}
				while($len_ori=~s/^([^\,]+)\,//)
				{
					my $len_=$1;
					if ($len=~/\w/) { $len.=","; }
					$len.="$len_";
				}
				print OUT qq!$out\t$len\t$pos\n!;
			}
		}
	}
	close(OUT);
	
	# [MANOR] Hardcoding done here:
	open (OUT,">proteome.bed");
	
	open (LOG,">$filename.bed.log");
	if (open (IN,"$filename"))
	{
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^LOCUS/) 
			{ 
				if ($line=~/^LOCUS\s*([^\s]+)/) 
				{ 
					my $protein=$1; 
					my $gene="";
					my $mrna_name="";
					my $mrna_protein_start="";
					my $mrna_protein_end="";
					my $definition_cont=0;
					while ($line and $line!~/^\/\//)
					{
						if ($line=~/^\s*\/coded_by=\"([^\"]+)\"$/) 
						{
							my $line_=$1;
							if ($line_=~s/^\s*complement\((.*)\)\s*$/$1/) { print  LOG qq!Warning: $line\n!;}
							if ($line_=~s/\<//g) { print  LOG qq!Warning: Removed '<': $line\n!;}
							if ($line_=~s/\>//g) { print  LOG qq!Warning: Removed '>': $line\n!;}
							if ($line_=~/^([^\:]+):([0-9]+)\.\.([0-9]+)$/) 
							{ 
								$mrna_name=$1; 
								$mrna_protein_start=$2; 
								$mrna_protein_end=$3; 
								$mrna_name=~s/\..*$//;
								$protein_transcript{$protein}=$mrna_name;
							}
							else
							{
								print LOG qq!Error parsing coded_by: $line\n!;
							}
						}
						if ($line=~/^\s*\/gene=\"([^\"]+)\"$/) 
						{ 
							$gene=$1; 
							$genes{$protein}=$gene;
						}
						if ($definition_cont==1)
						{
							if ($line=~/^\s*(.*)$/)
							{
								$description{$protein}.=" $1";
								if ($description{$protein}=~/\.\s*$/)
								{
									$definition_cont=0;
								}
							}
						}
						if ($line=~/^DEFINITION\s*(.*)$/)
						{
							$description{$protein}=$1;
							if ($description{$protein}!~/\.\s*$/)
							{
								$definition_cont=1;
							}
						}
						if ($line=~/^VERSION\s*(.*)$/)
						{
							my $gi=$1;
							if ($gi=~/GI:([0-9]+)/) 
							{
								$gi=$1;
								$gis{$protein}=$gi;
							}
						}
						$line=<IN>;
						chomp($line);
					}
					if ($mrna_name=~/\w/ and $mrna_protein_start=~/\w/ and $mrna_protein_end=~/\w/) 
					{
						if ($alignment_chr{$mrna_name}=~/\w/)
						{
							if ($alignment_strand{$mrna_name}=~/\-/)
							{
								my $len_sum=0;
								my $temp=$alignment_seg_len{$mrna_name};
								while($temp=~s/^([^\,]+)\,//)
								{
									my $len=$1;
									$len_sum+=$len;
								}
								my $temp_=$mrna_protein_start;
								$mrna_protein_start=$len_sum-$mrna_protein_end+1;
								$mrna_protein_end=$len_sum-$temp_+1;
							}
							my $started=0;
							my $last=0;
							my $len_sum=0;
							my @exon_len=();
							my @exon_pos=();
							my $exon_num=0;
							my $min=0;
							my $max=0;
							my $temp=$alignment_seg_len{$mrna_name};
							my $temp_=$alignment_seg_pos{$mrna_name};
							while($temp=~s/^([^\,]+)\,//)
							{
								my $len=$1;
								if ($last==0)
								{
									if ($temp_=~s/^([^\,]+)\,//)
									{
										my $pos=$1;
										$len_sum+=$len;
										#print LOG qq!$alignment_chr{$mrna_name}\t$min\t$max\t$gene-$mrna_name-$protein\t$len_sum\t$mrna_protein_start\t$mrna_protein_end\n!;
										if ($mrna_protein_start<=$len_sum and $started==0)
										{
											$started=1;
											$pos+=$len-($len_sum-$mrna_protein_start)-1;
											$min=$pos;
											$len=$len_sum-$mrna_protein_start+1;
										}
										if ($mrna_protein_end<=$len_sum and $started==1)
										{
											$last=1;
											$len=$mrna_protein_end-($len_sum-$len);
											$max=$pos+$len;
										}
										if ($started==1)
										{
											if ($exon_num==0 or ($pos-$min)>0)
											{
												$exon_len[$exon_num]=$len;
												$exon_pos[$exon_num]=$pos-$min;
												$exon_num++;
												$max=$pos+$len;
											}
										}
									} else { print LOG qq!Error parsing position\n! }
								}
							}
							print OUT qq!$alignment_chr{$mrna_name}\t$min\t$max\t$protein\t1000\t$alignment_strand{$mrna_name}\t$min\t$max\t0\t$exon_num!;
							print OUT qq!\t!;
							for(my $i=0; $i<$exon_num;$i++) { if ($i>0) { print OUT ","; } print OUT "$exon_len[$i]"; }
							print OUT qq!\t!;
							for(my $i=0; $i<$exon_num;$i++) { if ($i>0) { print OUT ","; } print OUT "$exon_pos[$i]"; }
							print OUT qq!\n!;
						} else { print LOG qq!Error finding alinment: $protein $mrna_name\n! }
					} else { print LOG qq!Error finding mRNA: $protein\n! }
				} else { print LOG qq!Error finding protein name: $line\n! }
			}
			else
			{
				$line=<IN>;
			}
		}
		close(IN);	
	}
	if (open (OUT_DESC,">$filename-descriptions.txt"))
	{
		foreach my $protein (sort keys %description)
		{
			print OUT_DESC qq!$protein\t$description{$protein}\n!;
		}
		close(OUT_DESC);
	}
	if (open (OUT_GENE,">$filename-genes.txt"))
	{
		foreach my $protein (sort keys %genes)
		{
			print OUT_GENE qq!$protein\t$gis{$protein}\t\t$genes{$protein}\n!;
		}
		close(OUT_GENE);
	}
	if (open (OUT,">$filename-protein-transcript.txt"))
	{
		foreach my $protein (sort keys %protein_transcript)
		{
			print OUT qq!$protein\t$protein_transcript{$protein}\n!;
		}
		close(OUT);
	}
}