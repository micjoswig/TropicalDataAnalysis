# ./computations.rules

use application "tropical";
use Benchmark qw(:all);

# This procedure returns a random set of cardinality {$n_samples} using numbers from 0 to {$upper_bound} - 1.
sub random_indices {
    my ($upper_bound, $n_samples) = @_;
    my $random_indices = new Set();
    while ($random_indices -> size() < $n_samples) {
        my $random_number = int(rand($upper_bound));
        $random_indices -> collect($random_number);
    }
    return $random_indices;
}

# The three parameters are {$input_filename}, {$max_n_trees}, and {$n_iterations} representing the name of the file containing 300 randomly generated trees,
# the upper bound of the number of trees in a computation, and the number of iterations performed for a computation, respectively.
# For every integer {$i} from 1 to {$max_n_trees}, the procedure computes the times needed to compute an asymmetric tropical Fermat--Weber set of {$i} trees from {$input_filename}.
# The timings respresent averages of {$n_iterations} iterations. In each iteration, {$i} trees are randomly sampled from {$input_filename}.
sub average_exact_timings_asym {
    my ($input_filename, $max_n_trees, $n_iterations) = @_;
    my $arr = graph::read_trees($input_filename);
    my $equi_arr = makeEquidistant($arr);
    my $M = vectorize_trees($equi_arr);
    my @timings_array = ();
    
    for (my $i = 1; $i <= $max_n_trees; ++$i) {
        my $total_time = 0.0;
        for (my $j = 0; $j < $n_iterations; ++$j) {
            my $indices = random_indices(300, $i);
            my $V = new Matrix($M -> minor($indices, All));
            
            my $t0 = Benchmark -> new;
            fw_set($V);
            my $t1 = Benchmark -> new;
            
            my $timing = timediff($t1, $t0);
            my $string_timing = timestr($timing);
            if ($string_timing =~ /(\d+.\d\d) CPU\)$/) {
                my $iteration_time = new Float($1);
                $total_time += $iteration_time;
            }
        }
        push(@timings_array, $total_time / $n_iterations);
    }
    return @timings_array;
}

# The following procedure prints in the file specified by {$output_filename} a code for a LaTeX tabular which
# contains the timings in seconds for computing asymmetric tropical Fermat---Weber sets for up to {$max_n_trees} trees
# having between {$min_n_leaves} and {$max_n_leaves} leaves.
# The timings are averages over {$n_iterations} iterations.
# parameters ($max_n_trees, $min_n_leaves, $max_n_leaves, $n_iterations, $output_filename)
sub table_exact_asym {
    my ($max_n_trees, $min_n_leaves, $max_n_leaves, $n_iterations, $output_filename) = @_;
    open(FH, '>', $output_filename) or die $!;
    
    print FH "\\begin{tabular}{c";
    for (my $i = 1; $i <= $max_n_trees; ++$i) {
        print FH "r";
    } 
    print FH "}\n \\toprule\n Leaves\$\\backslash\$Trees ";
    for (my $i = 1; $i <= $max_n_trees; ++$i) {
        print FH " \& \\multicolumn{1}{c}{", $i, "}";
    }
    print FH "\\\\\n \\midrule\n";
    
    for (my $i = $min_n_leaves; $i <= $max_n_leaves; ++$i) {
        my $input_filename = './300trees-' . $i .'leaves.txt';
        my @timings_array = average_exact_timings_asym($input_filename, $max_n_trees, $n_iterations);
        print FH " ", $i, " \& ";
        print FH join(" \& ", @timings_array);
        print FH " \\\\ \n";
    }
 
    print FH "\\bottomrule \n\\end{tabular} \n";
    close(FH);
}

# The three parameters are {$input_filename}, {$max_n_trees}, and {$n_iterations} representing the name of the file containing 300 randomly generated trees,
# the upper bound of the number of trees in a computation, and the number of iterations performed for a computation, respectively.
# For every integer {$i} from 1 to {$max_n_trees}, the procedure computes the times needed to compute an asymmetric tropical Fermat--Weber set of {$i} trees from {$input_filename}.
# The timings respresent averages of {$n_iterations} iterations. In each iteration, {$i} trees are randomly sampled from {$input_filename}.
sub average_exact_timings_sym {
    my ($input_filename, $max_n_trees, $n_iterations) = @_;
    my $arr = graph::read_trees($input_filename);
    my $equi_arr = makeEquidistant($arr);
    my $M = vectorize_trees($equi_arr);
    my @timings_array = ();
    
    for (my $i = 1; $i <= $max_n_trees; ++$i) {
        my $total_time = 0.0;
        for (my $j = 0; $j < $n_iterations; ++$j) {
            my $indices = random_indices(300, $i);
            my $V = new Matrix($M -> minor($indices, All));
            
            my $t0 = Benchmark -> new;
            sym_fw_set($V);
            my $t1 = Benchmark -> new;
            
            my $timing = timediff($t1, $t0);
            my $string_timing = timestr($timing);
            if ($string_timing =~ /(\d+.\d\d) CPU\)$/) {
                my $iteration_time = new Float($1);
                $total_time += $iteration_time;
            }
        }
        push(@timings_array, $total_time / $n_iterations);
    }
    return @timings_array;
}


# The following procedure prints in the file specified by {$output_filename} a code for a LaTeX tabular which
# contains the timings in seconds for computing symmetric tropical Fermat---Weber sets for up to {$max_n_trees} trees
# having between {$min_n_leaves} and {$max_n_leaves} leaves.
# The timings are averages over {$n_iterations} iterations.
# parameters ($max_n_trees, $min_n_leaves, $max_n_leaves, $n_iterations, $output_filename)
sub table_exact_sym {
    my ($max_n_trees, $min_n_leaves, $max_n_leaves, $n_iterations, $output_filename) = @_;
    open(FH, '>', $output_filename) or die $!;
    
    print FH "\\begin{tabular}{c";
    for (my $i = 1; $i <= $max_n_trees; ++$i) {
        print FH "r";
    } 
    print FH "}\n \\toprule\n Leaves\$\\backslash\$Trees ";
    for (my $i = 1; $i <= $max_n_trees; ++$i) {
        print FH " \& \\multicolumn{1}{c}{", $i, "}";
    }
    print FH "\\\\\n \\midrule\n";
    
    for (my $i = $min_n_leaves; $i <= $max_n_leaves; ++$i) {
        my $input_filename = './300trees-' . $i .'leaves.txt';
        my @timings_array = average_exact_timings_sym($input_filename, $max_n_trees, $n_iterations);
        print FH " ", $i, " \& ";
        print FH join(" \& ", @timings_array);
        print FH " \\\\ \n";
    }
 
    print FH "\\bottomrule \n\\end{tabular} \n";
    close(FH);
}

# The three parameters are {$input_filename}, {$max_n_trees}, and {$n_iterations} representing the name of the file containing 300 randomly generated trees,
# the upper bound of the number of trees in a computation, and the number of iterations performed for a computation, respectively.
# For every integer {$i} from 1 to {$max_n_trees}, the procedure computes the times needed to compute an asymmetric tropical Fermat--Weber set of {$i} trees from {$input_filename}.
# The timings respresent averages of {$n_iterations} iterations. In each iteration, {$i} trees are randomly sampled from {$input_filename}.
# Note that this code calls mcf only if it is installed.
sub average_timings_mcf {
    my ($input_filename, $max_n_trees, $n_iterations) = @_;
    my $arr = graph::read_trees<Float>($input_filename);
    my $equi_arr = makeEquidistant($arr);
    my $M = vectorize_trees($equi_arr);
    my @timings_array = ();
    
    for (my $i = 50; $i <= $max_n_trees; $i += 50) {
        my $total_time = 0.0;
        for (my $j = 0; $j < $n_iterations; ++$j) {
            # for $i = 300, the random sampling is useless
            if ($i == 300) {
                my $t0 = Benchmark -> new;
                fw_set($M);
                my $t1 = Benchmark -> new;
                
                my $timing = timediff($t1, $t0);
                my $string_timing = timestr($timing);
                if ($string_timing =~ /(\d+.\d\d) CPU\)$/) {
                    my $iteration_time = new Float($1);
                    $total_time += $iteration_time;
                }
                next;
            }
            
            my $indices = random_indices(300, $i);
            my $V = new Matrix<Float>($M -> minor($indices, All));
            
            my $t0 = Benchmark -> new;
            fw_set($V);
            my $t1 = Benchmark -> new;
            
            my $timing = timediff($t1, $t0);
            my $string_timing = timestr($timing);
            if ($string_timing =~ /(\d+.\d\d) CPU\)$/) {
                my $iteration_time = new Float($1);
                $total_time += $iteration_time;
            }
        }
        push(@timings_array, $total_time / $n_iterations);
    }
    return @timings_array;
}

# The following procedure prints in the file specified by {$output_filename} a code for a LaTeX tabular which
# contains the timings in seconds for computing symmetric tropical Fermat---Weber sets for up to {$max_n_trees} trees
# having between {$min_n_leaves} and {$max_n_leaves} leaves.
# The timings are averages over {$n_iterations} iterations.
# Note that this code calls mcf only if it is installed.
# parameters ($max_n_trees, $min_n_leaves, $max_n_leaves, $n_iterations, $output_filename)
sub table_mcf {
    my ($max_n_trees, $min_n_leaves, $max_n_leaves, $n_iterations, $output_filename) = @_;
    open(FH, '>', $output_filename) or die $!;
    
    print FH "\\begin{tabular}{c";
    for (my $i = 50; $i <= $max_n_trees; $i += 50) {
        print FH "r";
    } 
    print FH "}\n \\toprule\n Leaves\$\\backslash\$Trees ";
    for (my $i = 50; $i <= $max_n_trees; $i += 50) {
        print FH " \& \\multicolumn{1}{c}{", $i, "}";
    }
    print FH "\\\\\n \\midrule\n";
    
    for (my $i = $min_n_leaves; $i <= $max_n_leaves; $i += 5) {
        my $input_filename = './300trees-' . $i .'leaves.txt';
        my @timings_array = average_timings_mcf($input_filename, $max_n_trees, $n_iterations);
        print FH " ", $i, " \& ";
        print FH join(" \& ", @timings_array);
        print FH " \\\\ \n";
    }
 
    print FH "\\bottomrule \n\\end{tabular} \n";
    close(FH);
}

# The following procedure reads trees from a given input file, given by {$input_filename}, and computes
# timings of computing asymmetric tropical Fermat--Weber sets using mcf. The timings are obtained as
# averages over {$n_iterations}. The average times are printed each on a line in the file with the path
# given by {$output_filename}. 
# Note that this code calls mcf only if it is installed.
# parameters ($input_filename, $output_filename, $min_n_trees, $max_n_trees, $n_iterations)
sub timings_mcf {
    my ($input_filename, $output_filename, $min_n_trees, $max_n_trees, $n_iterations) = @_;
    open(FH, '>', $output_filename) or die $!;
    
    my $arr = graph::read_trees<Float>($input_filename);
    my $equi_arr = makeEquidistant($arr);
    my $M = vectorize_trees($equi_arr);
    
    for (my $i = $min_n_trees; $i <= $max_n_trees; ++$i) {
        my $total_time = 0.0;
        print FH $i, " ";
        for (my $j = 0; $j < $n_iterations; ++$j) {
            my $indices = random_indices(300, $i);
            my $V = new Matrix<Float>($M -> minor($indices, All));
            
            my $t0 = Benchmark -> new;
            fw_set($V);
            my $t1 = Benchmark -> new;
            
            my $timing = timediff($t1, $t0);
            my $string_timing = timestr($timing);
            if ($string_timing =~ /(\d+.\d\d) CPU\)$/) {
                my $iteration_time = new Float($1);
                $total_time += $iteration_time;
            }
        }
        print FH ($total_time / $n_iterations);
        print FH "\n";
    }
    close(FH);
}
