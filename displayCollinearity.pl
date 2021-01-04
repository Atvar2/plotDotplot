#! /usr/bin/perl -w

=head1 Name

 displayCollinearity.pl 

=head1 Description

 This program can show the collinearity relationship between two species with dotplot,
 basing on identified blocks by mscan or alignment of two genome by mummer.

=head1 Version

 Author: chenjhbio, chenjunhui@genomics.cn
 Version: 1.0
 Date: 20200428
 Note:
 1.if set --identity=0(deafault),plot the "+" and "-" direction synteny dots with red and blue color,
 otherwise,it will show aligment between two genomes by a range color.
 2. --reorder option,which imply your aliment may be reorder in order to obtain a better display.

=head1 Usage

 perl displayCollinearity.pl  <syn_chr_file> <chr1_file> <chr2_file> <svg_out> [options]
 
 <blockFile>:   The synteny relationshiop or aligment between two genome. 
 <queryLength>:  The species chromosomes length file plot in the X axis.
 <referenceLength>:The species chromosomes length file plot in the Y axis.
 <svg_out>:      The output svg file.
 --identity      Add identity to the figure if display the alignment.  [default: 0]
 --spec1 <str>:  Set the X axis title. [default:species1]
 --spec2 <str>:  Set the Y axis title. [default:species2]
 --r     <int>:  Set the radius of each syntenic dot. [default:3]
 --unit  <int>:  Set the scale on the X/Y axis. [default:20]
 --rx    <int>:  Set the resolution of X axis. [default:200000]
 --ry    <int>:  Set the resolution of Y axis. [default:200000]
 --reorder:      reorder chromosome of x axis basing on chromosome on y axis.
 --orderaxis:    axis needed to be ordered, default: x.
 --png:          Convert the svg file into png format.

=head1 Example

 perl displayCollinearity.pl  blockFile  [queryLength]  [referenceLength]  output.svg --spec1 specie --spec2  specie  --r 2 -- unit 20 --reorder orderFile --orderaxis x  --png


=cut

use strict;
use FindBin qw($Bin);
use lib "plotDotplot/lib";
use SVG;
use Font;
use Getopt::Long;
use Data::Dumper;
my ($spec1,$spec2,$iden,$num,$maxiden,$r,$png,$Axis_unit,$reorder,$axis,$resol_x,$resol_y);
GetOptions(
		"spec1:s"=>\$spec1,
		"spec2:s"=>\$spec2,
		"r:s"=>\$r,
		"identity:i"=>\$iden,
		"num:i" =>\$num,
		"max:i" =>\$maxiden,
		"unit:i"=>\$Axis_unit,
		"png"=>\$png,
                "reorder:s"=>\$reorder,
		"orderaxis:s"=>\$axis,
		"rx:i"=>\$resol_x,
		"ry:i"=>\$resol_y,
);
$spec1 ||= 'species1';
$spec2 ||= 'species2';
$r ||= 3;
$iden ||=0;
$maxiden ||=100;
$num     ||=10;
$axis ||="y";
my $start ||=85;
$Axis_unit ||= 20;
$resol_x ||= 200000;
$resol_y ||= 200000;
my $mark ||=3;
die `pod2text $0` if (@ARGV == 0);

my $syn_chr_file = shift;
my $chr1_file = shift;
my $chr2_file = shift;
my $colorfile =shift;
my $svg_out_file = shift;	
my @colors = &getColor($colorfile);
my @Markcolor   =&getMarkColor(\@colors,$mark);
#my $svg2xxx ="/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_evolution/Tools/svg2xxx_release/svg2xxx";
my %stepnum=();my %step=();my %match=();my @pos;
my %step255=();
my %stepround=();my @gradualColors;
if($iden ==1){
	&dividcolor($maxiden,$num);
	@gradualColors=&Producegradual("255,0,0","204,255,204")
}

my (%hash,%dire);		

open SYN,"$syn_chr_file" || die "$!";
while (<SYN>)
{
	my @temp = split /\s+/,$_;
	my ($s1,$start1,$end1,$s2,$start2,$end2,$direct) = @temp[0,2,3,4,6,7,8];
	($start1,$end1) = ($end1,$start1) if ($start1 > $end1);
	($start2,$end2) = ($end2,$start2) if ($start2 > $end2);
	push @{$hash{$s1}{$s2}},[$start1,$start2,$direct];
}
close SYN;

############  chr1 x axis  information
open C1,"$chr1_file" || die "$!";
my $tot_len1;
my %len1;
while (<C1>){
	chomp;
	my ($name,$len) = split /\s+/,$_;
	$name = substr($name,0,20);
	$name =~ s/\D+//g;###########################
#	$name = scalar($name);
	$len1{$name} = $len;
	$tot_len1 += $len;
}
close C1;

############# chr2 ix axis information
open C2,"$chr2_file" || die "$!";
my $tot_len2;
my %len2;
while (<C2>){
	chomp;
	my ($name,$len) = split /\s+/,$_;
	$name = substr($name,0,20);
	$name =~ s/\D+//g;#######################
	$len2{$name} = $len;
	$tot_len2 += $len;
}
close C2;
############################### deal with chr order basing on alignment
my @orderChr;
if($reorder){
        &getMatchpos(\$reorder);
	if($axis eq 'x'){
		@orderChr=keys %len1;
	}
	else{
		@orderChr=keys %len2;
	}

}
else{
	        @pos=sort {$a<=>$b} keys %len1;
		@orderChr=@pos;
}

if(scalar(@orderChr) !=scalar(@pos)){
	print "Number of chromosomes in order file are not consistent with the orginal, please check it!!!";
	exit(1);
}
##############################    control area
#my ($resol_x,$resol_y)=(200000,200000);

my ($gen_width,$gen_height) = ($tot_len1/$resol_x,$tot_len2/$resol_y);
my ($x_blank,$y_blank) = ($gen_width*0.12,$gen_height*0.15);
my ($width,$height) = ($gen_width + 2*$x_blank,$gen_height + 2*$y_blank);
my ($t_x,$t_y) = ($x_blank*1.2,$y_blank*0.6);
#my ($t_x,$t_y) = ($x_blank,$y_blank);
my ($x_text_v,$y_text_v) = ($x_blank*0.5,$y_blank*0.8);
my $stroke = 'black';
my $stroke_width = $gen_width/400;
my ($font_family,$font_fam2) = ('Courier','Courier');
my $font_size = $gen_width/55;
### draw svg 
my $svg = SVG->new(width=>$width,height=>$height);
$svg ->rect('id'=>'bg_groud','x'=>0,'y'=>0,width=>$width,height=>$height,'stroke'=>'white','fill'=>'white');
my $L1 = $svg->group('stroke'=>$stroke,'fill'=>'none','stroke-width'=>$stroke_width);
my $L2 = $svg->group('stroke'=>$stroke,'stroke-width'=>$stroke_width*0.5);

my $T1 = $svg->group('text-anchor'=>'middle','font-family'=>$font_family,'font-size'=>$font_size,'stroke'=>$stroke,'fill'=>$stroke);
my $T2 = $svg->group('text-anchor'=>'middle','font-family'=>$font_fam2,'font-size'=>$font_size*0.8,'stroke'=>$stroke,'fill'=>$stroke);
my $T3 = $svg->group('text-anchor'=>'middle','font-family'=>$font_family,'font-size'=>$font_size*1.2,'stroke'=>'black','fill'=>'black','font-style'=>'italic');
$L1 -> rect(id=>'out','x'=>$t_x,'y'=>$t_y,'width'=>$gen_width,'height'=>$gen_height);
#$L2 -> line('x1'=>$t_x,'y1'=>$t_y+$gen_height,'x2'=>$t_x+$gen_width,'y2'=>$t_y)if(defined $mirror);
##########################  spec name
my $spec2_text_x = $t_x - $x_blank*0.9;
my $spec2_text_y = $t_y + $gen_height/2;
$T3 -> text('id','spec2','x',$spec2_text_x,'y',$spec2_text_y,'transform',"rotate(-90,$spec2_text_x,$spec2_text_y)")->cdata($spec2);

my $spec1_text_x = $t_x+$gen_width/2;
my $spec1_text_y = $t_y+$gen_height+$y_blank*0.5;   # x axis  1.12 -> 0.5
$T3 -> text('id','spec1','x',$spec1_text_x,'y',$spec1_text_y)->cdata($spec1);

##########################   vertical line
my $x_add = 0;
my %x_p_hash; #####  store the x point of x_chr
#print Dumper(%len1);
foreach my $chrpos(@pos){
	my $text_temp = $x_add + $t_x;
	$x_p_hash{$chrpos} = $text_temp;
	$x_add += $len1{$chrpos}/$resol_x;
	my $chr_x = $x_add + $t_x;
	my $y1 = $t_y + $gen_height;
	my $y2 = $t_y;
	$L2->line('id'=>"x_$chrpos",'x1'=>$chr_x,'y1'=>$y1,'x2'=>$chr_x,'y2'=>$y2);

	my $text_x = ($text_temp + $chr_x)/2;
	my $text_y = $y1 + $y_text_v*0.35;  # 0.5 ->0.3
	$T1->text('id'=>"x_text_$chrpos",'x'=>$text_x,'y'=>$text_y,'transform'=>"rotate(-30,$text_x,$text_y)")->cdata("chr$chrpos");
}

###########################   cross line
my $y_add = 0;
my %y_p_hash;  #####  store the y point of y_chr
foreach (sort {$a<=>$b} keys %len2){
	my $text_temp = $t_y + $gen_height - $y_add;
	$y_p_hash{$_} = $text_temp;
#print "$_\tuuuuuuu\n";
	$y_add += $len2{$_}/$resol_y;
	my $chr_y = $t_y +$gen_height - $y_add;
	my $x2 = $t_x + $gen_width;
	$L2->line('id'=>"y_$_",'x1'=>$t_x,'y1'=>$chr_y,'x2'=>$x2,'y2'=>$chr_y);

	my $text_x = $t_x - $x_text_v*0.9; #1->0.9
	my $text_y = ($text_temp + $chr_y)/2;
	$T1->text('id'=>"y_text_$_",'x'=>$text_x,'y'=>$text_y,'transform'=>"rotate(-30,$text_x,$text_y)")->cdata("chr$_");
}

############################  dot plot
my $circle_r = $svg->group('id'=>'circle_r','stroke'=>'red','fill'=>'red');
my $circle_b = $svg->group('id'=>'circle_b','stroke'=>'blue','fill'=>'blue');
my $circle   =$svg->group('id'=>'circle');
foreach my $s1(keys %hash){
	my $s1_key; if ($s1 =~ /(\d+)/) { $s1_key = $1; }
	foreach my $s2(keys %{$hash{$s1}}){
		my $s2_key; if ($s2 =~ /(\d+)/) { $s2_key = $1; }
		for(my $i=0;$i<@{$hash{$s1}{$s2}};$i++)
        {
			my $sta1 = ${$hash{$s1}{$s2}}[$i][0]; ### y ver
			my $sta2 = ${$hash{$s1}{$s2}}[$i][1]; ### x ver
		        my $dire = ${$hash{$s1}{$s2}}[$i][2]; ### x ver
			my $p_y=0;my $p_x=0;
			if (!defined $y_p_hash{$s1_key} or !defined $x_p_hash{$s2_key}){
				print "warnings:chromoses in xaxis or yaixs may not plotted for alignment:$s1\t$s2\t$sta1, please check it!!!\n";
			}
			else{
				$p_y = $y_p_hash{$s1_key} - $sta1/$resol_y;
				$p_x = $x_p_hash{$s2_key} + $sta2/$resol_x;
			}
			if ($iden ==0){
				$circle_r->circle('cx'=>$p_x,'cy'=>$p_y,'r'=>$r)if($dire eq 'plus');
				$circle_b->circle('cx'=>$p_x,'cy'=>$p_y,'r'=>$r)if($dire eq 'minus');
			}
			elsif($iden ==1){
				my $color=&num2color($dire);
				if (! $color){
					print "$dire\n";
				}
				if(!($color)){$color='red';}
				$color="rgb".$color;
				$circle->circle('cx'=>$p_x,'cy'=>$p_y,'r'=>$r,'fill'=>$color,'stroke'=>$color);
			}
		}
	}
}

########################## x mark
#my $Axis_unit = 20;
my $gap = 1000*1000*$Axis_unit/$resol_x;
my $m_t_x = $t_x;
my $flag = 0;
my $num_x = 0;
for (my $i=$gap;$i<$gen_width;$i+=$gap){
	$m_t_x += $gap;
	my $y1 = $t_y + $gen_height;
	my $y2 = $y1 + $y_blank*0.1;     #0.12->0.1  axis mark  change 20190430
	$L1->line('id'=>"x_mark_$i",'x1'=>$m_t_x,'y1'=>$y1,'x2'=>$m_t_x,'y2'=>$y2);

	if ($flag % 2){
#		my $text_x = $m_t_x;
		$num_x += 2*$Axis_unit;
		my $text_y = $y2 + $y_blank * 0.26;
		#$T2->text('id'=>"x_m_t_$i",'x'=>$m_t_x,'y'=>$text_y)->cdata("$num_x".'Mb');
	}
	$flag += 1;
}

######################## y mark
my $m_t_y = $t_y + $gen_height;
$flag = 0;
my $num_y = 0;
for (my $i=$gap;$i<$gen_height;$i+=$gap){
	$m_t_y -= $gap*$resol_x/$resol_y;
	my $x1 = $t_x;
	my $x2 = $t_x - $x_blank*0.10;
	$L1->line('id'=>"y_mark_$i",'x1'=>$x1,'y1'=>$m_t_y,'x2'=>$x2,'y2'=>$m_t_y);
	if ($flag % 2){
		$num_y += 2*$Axis_unit;
		my $text_x = $x2 - $x_blank * 0.05;
		#$T2->text('id'=>"y_m_t_$i",'x'=>$text_x,'y'=>$m_t_y,'transform'=>"rotate(-90,$text_x,$m_t_y)")->cdata("$num_y".'Mb');
	}
	$flag += 1;
}
#################  color code
#@Markcolor
my $x_code=$gen_width+$x_blank*1.4;my $y_code=$y_blank*1.2;
my $Cwidth=40;my $Cheight=sprintf("%d",$gen_height/4);
my $eachheight=$Cheight/4;
my $lg2 = $svg->gradient(-type=>'linear',
                                id => 'code',
                                x1=>"0%",
                                y1=>"0%",
                                x2=>"0%",
                                y2=>"100%",
                                spreadMethod=>"pad",
                                gradientUnits=>"userSpaceOnUse");

        $lg2->stop(offset=>"0%",
            style=>{'stop-color'=>"red",'stop-opacity'=>"1"});
	 $lg2->stop(offset=>"25%",
            style=>{'stop-color'=>"rgb$colors[4]",'stop-opacity'=>"1"});
        $lg2->stop(offset=>"50%",
            style=>{'stop-color'=>"rgb$colors[8]",'stop-opacity'=>"1"});
	 $lg2->stop(offset=>"75%",
            style=>{'stop-color'=>"rgb$colors[12]",'stop-opacity'=>"1"});
        $lg2->stop(offset=>"100%",
            style=>{'stop-color'=>"green",'stop-opacity'=>"1"});

$svg->rectangle(
    x => $x_code,
    y => $y_code,
    rx=> 20,
    ry=> 20,
    width  => $Cwidth,
    height => $Cheight,
    style => {
        'fill'           => 'url(#code)',
        'stroke'         => 'url(#code)',
        'stroke-width'   =>  0.1,
        'stroke-opacity' =>  1,
        'fill-opacity'   =>  1,
    },
);
############# code marks
my $x_mark=$gen_width+$x_blank*1.4+40;my $y_mark=$y_blank*1.2;
my $x_text=$gen_width+$x_blank*1.4+100;my $y_text=$y_blank*1.2;
my $tmpy_mark=$y_mark;my $tmpy_text=$y_text;
my $T4 = $svg->group('text-anchor'=>'middle','font-size'=>40,'stroke'=>'black','fill'=>'black');
my $L6 = $svg->group('stroke' =>'black','stoke-width'=>1,'stroke-opacity' => 1,'fill-opacity'=>1);
for (my $i=0;$i<5;$i++){
	my $text=$i*25;
	$y_mark=$tmpy_mark+$eachheight*$i;
	$y_text=$tmpy_text+$eachheight*$i;
	$L6->rectangle('id' =>"rect_$i",'x'=>$x_mark,'y'=>$y_mark,width=>10,height=>2);
	$T4->text('id' =>"marker_$i",'x'=>$x_text,'y'=>$y_text)->cdata("$text%"); 
	
}

#############get color #######
sub dividcolor(){
	my $max=shift;
	my $start=85;
	$max ||=100;my $step=$max/$num;
	my $step255=255/$num;
	for(0..$num){
		$stepnum{$_}=$step*$_;
		$step{$_}=$step255*$_;
		$stepround{$_}=int($step255*$_ + 0.5); # same as %step, but in pround format.
	}
}

sub  num2color(){
	my $number=shift;
	my $j=0;
	my $color="";
	for(my $i=$start;$i<=100;$i++){
		if($number >=$i && $number<=$i+1){
			$color=$colors[$j];
		}
		if($number == $start){$color=$colors[0];}
		if($number == 100){$color=$colors[-1];}
		$j++;
	}
	return $color;
}

sub getMarkColor(){
	my $p=shift;
	my $m=shift;
	my @mark;
	my $len=scalar(@{$p});
	for(my $i=0;$i <=$len;$i++){
		if($i % $m == 0){
			push @mark,$p->[$i];
		}
	}
	return  @mark;
}

sub getColor(){
	my $file=shift;
	my @array;
	open IN, $file;
        my $line=<IN>;chomp($line);
	while($line=~/\"(\S+)\"/g){
		my $color=$1;
	        $color=~s/#//g;
		my $r=hex(substr($color,0,2));
	        my $g=hex(substr($color,2,2));
                my $b=hex(substr($color,4,2));
		push  @array,"($r,$g,$b)";
	}
	return  @array;
}

sub Producegradual(){
	my $RGB1=shift;
	my $RGB2=shift;
	my $Tstep=shift;
	$Tstep   ||=10;
	my($RA,$GA,$BA)=split(",",$RGB1);
	my($RB,$GB,$BB)=split(",",$RGB2);
	my @gradualcolor;
	for (my $step=1;$step <= $Tstep;$step++){
		my $Rstep=$RA+($RB-$RA)*$Tstep/$step;
		my $Gstep=$GA+($GB-$GA)*$Tstep/$step;
		my $Bstep=$BA+($BB-$BA)*$Tstep/$step;
		push @gradualcolor,"$Rstep,$Gstep,$Bstep";
	}
	return  @gradualcolor;
}
##################### get match chr pos ######################
sub getMatchpos($){
	#get match chr pos if chromosomes do not match exactly on 
	#diagonal;
	my $p=shift;
	open IN,$$p;
        while(<IN>){
                chomp;
                my @array=split(/\t/);
                $array[0]=~s/\D+//g;
		$array[1]=~s/\D+//g;
                $match{$array[0]}=$array[1];
        }
        my @refChr=sort {$a <=> $b} keys(%match);
        for(@refChr){
                push @pos,$match{$_};
        }
	return 0;
}

##########################  print out
open OUT,">$svg_out_file" || die "$!";
print OUT $svg->xmlify();
close OUT;
