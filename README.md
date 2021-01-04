# plotDotplot
The kit was used for displaying dotplot of genome alignment basing on mscan and mummer
===
Copyright (C) 2014 Arthur Zwaenepoel
   Apache License Version 2.0, January 2004 http://www.apache.org/licenses/
  
 Beijing Genomics Institute, Animal and Plant research https://www.genomics.cn/

#The Kit was used for displaying genomic collinearity basing on result of [Mummer](http://mummer.sourceforge.net/)  or [Mscan](https://github.com/tanghaibao/mcscan). 

Installation
===
---
      The Kit was writted mainly by perl, and tested  with perl v5.10.1 on CentOS release 6.9.
      Install it, just simply execute as following.Then you can run it directly.
      ```git clone  https://github.com/chenjunhui/plotDotplot/edit/plotDotplot```
      
      Notably, the kit depend on SVG package of perl, thus, you should install it and `export PERL5LIB="/user/packages:$PERL5LIB"` before run it.
Usage
===
----
```
perl ../displayCollinearity.pl test0.syn  Actinidia_chinensis.Red5.len   kiwifruit.chr.len   colors.txt   test.svg --spec1 specie --spec2 specie --r 2 --unit 20  --identity 1
```

Parameters
===
---
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
