#!/usr/bin/perl -w

# prog1.1 
# Bruno Contreras-Moreira / Raul Mejia Pedroza / Raul Ossio Vela / Alberto Vazquez Salazar
# Nearest Neighbor dG calculator

use strict;
use Term::ANSIColor;


# global variables
my $T           = 37; # temperature(C)
my $windowL     = 15;  # window length, http://www.biomedcentral.com/1471-2105/6/1
my %NNparams    = ( 
	# SantaLucia J (1998) PNAS 95(4): 1460-1465.
	# [NaCl] 1M, 37C & pH=7 
	# H(enthalpy): kcal/mol	, S(entropy): cal/k�mol
	# stacking dinucleotides
	'AA/TT' , {'H',-7.9, 'S',-22.2},
	'AT/TA' , {'H',-7.2, 'S',-20.4},
	'TA/AT' , {'H',-7.2, 'S',-21.3},
	'CA/GT' , {'H',-8.5, 'S',-22.7},
	'GT/CA' , {'H',-8.4, 'S',-22.4},
	'CT/GA' , {'H',-7.8, 'S',-21.0},
	'GA/CT' , {'H',-8.2, 'S',-22.2},
	'CG/GC' , {'H',-10.6,'S',-27.2},
	'GC/CG' , {'H',-9.8, 'S',-24.4},
	'GG/CC' , {'H',-8.0, 'S',-19.9},
	# initiation costs
	'G'     , {'H', 0.1, 'S',-2.8 },
	'A'     , {'H', 2.3, 'S',4.1  },
	# symmetry correction
	'sym'   , {'H',   0, 'S',-1.4 } );

# Declaracion de tablas para guardar datos
my @tabla_dG = ();
my @tabla_E1 = ();
my @tabla_E2 = ();
my @tabla_Dn = ();
my @tabla_seq_15 = ();
my @tabla_promotor = ();

# Declaracion e inicializacion de arreglo para contar numero de promotores de todas las secuencias en un archivo.
my $cc =1;
while ($cc < 451)
{
	$tabla_promotor[$cc] = 0;
	$cc++;
}
# Variables cutoff / pueden modificarse para cambiar especificidad en busqueda de promotores
my $cutoff_1 = (3);
my $cutoff_2 = (-18);

# Variables que definen el inicio y tamano de las ventanas para calcular E1 y E2
my $vent_E1_start = 0;
my $vent_E1 = 50;

my $vent_E2_start = 50;
my $vent_E2 = 50;

# Variables para calcular posicion de acuerdo al 0 y contar el numero de secuencias leidas por el programa
my $pos =0;
my $cont_seq = 1;


my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";

print "# parameters: Temperature=$T\C Window=$windowL\n\n";

open(SEQ, $infile) || die "# cannot open input $infile : $!\n";
while(<SEQ>)
{
	if(/^(b\d{4}) \\ ([ATGC]+)/)
	{
		my ($name,$seq) = ($1,$2); 
		# Variables para almacenar fragmentos de la secuencia que se enviaran a las subrutinas calc_dG, calc_E1 y calc_E2, asi como sus retrornos
		my $seq_15;
		my $seq_E1;
		my $seq_E2;
	
		my $n_dG = 0;
		my $E1 = 0;
		my $E2 = 0;



		printf("sequence %s (%d nts) $cont_seq\n",$name,length($seq));
		my $n = 1;
		my $cont = 0;

		while($cont < (length($seq)-$windowL))
		{
			$seq_15 = substr($seq, $n, $windowL);
			$n_dG = duplex_deltaG($seq_15,$T);
			#print "SEQ_$n $seq_15  $n_dG\n";
			$tabla_dG [$n] = $n_dG;
			$tabla_seq_15 [$n] = $seq_15;
			$cont++;
			$n++;
			$cont_seq++;
		}
		
		# Bucle para enviar secuencias a calc_E1
		$n = 1;
		$cont = 0;

		while($cont < (length($seq)-$vent_E1-1))
		{
			$seq_E1 = substr($seq, $n + $vent_E1_start, $vent_E1);
			$E1 = calc_E1($seq_E1,$tabla_dG [$n]);
			#print "SEQ_$n E1: $seq_E1 $E1\n";
			$tabla_E1 [$n] = $E1;
			$cont++;
			$n++;
		}
		# Bucle para enviar secuencias a calc_E2
		$n = 1;
		$cont =0;

		while($cont < (length($seq)-$vent_E2-1))
		{
			$seq_E2 = substr($seq, $n + $vent_E2_start, $vent_E2);
			$E2 = calc_E2($seq_E2, $tabla_dG [$n]);
			#print "SEQ_$n E2: $seq_E2 $E2\n";
			$tabla_E2 [$n] = $E2;
			$cont++;
			$n++;
		}
		# Bucle para enviar valores de E1 y E2 a calc_Dn
		$n = 1;
		$cont =0;

		while($cont < (length($seq)-$vent_E2-1))
		{
			my $Dn = abs(calc_Dn($tabla_E1 [$n],$tabla_E2 [$n]));
			#print "SEQ_$n Dn: $Dn\n";
			$tabla_Dn [$n] = $Dn;
			$cont++;
			$n++;
		}
		$n = 1;
		$cont = 0;
		my $promotor = " ";
		

		while($cont < (length($seq)-$vent_E2-1))
		{
			if(($tabla_Dn[$n] > $cutoff_1) && ($tabla_E1[$n] > $cutoff_2))
			{
				$promotor = 1;
			}
			else
			{
				$promotor = 0;
			}
			$pos = $n-401;
			print "Seq: $n  Pos: $pos\t $tabla_seq_15[$n]\t";
			printf ("%5.5f",$tabla_dG[$n]); 
			print "\t";
			printf ("%5.5f",$tabla_E1[$n]); 
			print "\t";
			printf ("%5.5f",$tabla_E2[$n]); 
			print " \t";
			printf ("%5.5f",$tabla_Dn[$n]); 
			print "   \t";

			#$tabla_E1[$n]\t $tabla_E2[$n]\t $tabla_Dn[$n]\t\t";

			if($promotor == 1)
			{
				print color('red');
				print " $promotor\t";
				$tabla_promotor[$n]++;
			} 
			else
			{
				print color('black');
				print " $promotor\t";
			}	
			print color ('black');
			
			print "$tabla_promotor[$n]\n";

			$cont++;
			$n++;
		}
		print "\n\n"; 
	}
	
}
close(SEQ);
# Crear archivo para poder hacer graficas de la frecuencia con la que se predicen los promotores
my $file = "temp.txt";

unless(open FILE, '>'.$file) {die "\nNo se pudo crear archivo $file para guardar datos\n";}
my $cc2 =1;
while ($cc2 < 451)
{
	my $pos2 = $cc2-401;
	print FILE "$pos2\t$tabla_promotor[$cc2]\n";
	$cc2++;
}
close FILE;

# calculate NN free energy of a DNA duplex , dG(t) = (1000*dH - t*dS) / 1000
# parameters: 1) DNA sequence string; 2) Celsius temperature
# returns; 1) free energy scalar
# uses global hash %NNparams
sub duplex_deltaG
{
   	my ($seq,$tCelsius) = @_; 
	
	my ($DNAstep,$nt,$dG,$total_dG) = ('','',0,0);
	my @sequence = split(//,uc($seq));
	my $tK = 273.15 + $tCelsius;
	
	sub complement{ $_[0] =~ tr/ATGC/TACG/; return $_[0] }
	
	# add dG for overlapping dinculeotides
	for(my $n=0;$n<$#sequence;$n++) 
	{
			$DNAstep = $sequence[$n].$sequence[$n+1].'/'.
				complement($sequence[$n].$sequence[$n+1]);
			
			if(!defined($NNparams{$DNAstep}))
			{
				$DNAstep = reverse($DNAstep);
			}
			
			$dG = ((1000*$NNparams{$DNAstep}{'H'})-
					($tK*$NNparams{$DNAstep}{'S'}))
					/ 1000 ;
			
			$total_dG += $dG; 
	}
	
	# add correction for helix initiation
	$nt = $sequence[0]; # first pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) } 
	$total_dG += ((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'}))
					/ 1000; 
	
	$nt = $sequence[$#sequence]; # last pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) }
	$total_dG += ((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'}))
					/ 1000;
	
	return $total_dG;
}

# Calcula valores de E1
# parameters: 1) Secuencia de ADN; 2) Valor de dG para una n calculado en base a una determinada ventana
# returns: Suma de todos los valores dG en la secuencia dividida entre el numero de elementos de la secuencia

sub calc_E1
{
	my ($E1_seq,$dG) = @_;
	my $E1_cont = 0;
	my $E1_calc = 0;

	while($E1_cont < length($E1_seq))
	{
		$E1_calc = $E1_calc + $dG;
		$E1_cont++;
	}

	return ($E1_calc/$vent_E1);
}

# Calcula valores de E2
# parameters: 1) Secuencia de ADN; 2) Valor de dG para una n calculado en base a una determinada ventana
# returns: 1) Suma de todos los valores dG en la secuencia dividida entre el numero de elementos de la secuencia

sub calc_E2
{
	my ($E2_seq,$dG) = @_;
	my $E2_cont = 0;
	my $E2_calc = 0;

	while ($E2_cont < length($E2_seq))
	{
		$E2_calc = $E2_calc + $dG;
		$E2_cont ++;
	}
	return ($E2_calc/$vent_E2);
}
# Realiza la resta E1-E2
# parameters: 1) Valor de E1 para una determinada n 2) Valor de E2 para una determinada n
# returns: resta E1-E2
sub calc_Dn
{
	my ($E1_D,$E2_D) = @_;
	my $D_calc = $E1_D - $E2_D;

	return $D_calc;
}


