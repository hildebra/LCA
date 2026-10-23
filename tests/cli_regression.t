#!/usr/bin/env perl
# End-to-end regression checks using temporary inputs and outputs only.
# Run with: LCA_BIN=./LCA prove tests/cli_regression.t
use strict;
use warnings;
use Cwd qw(abs_path);
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;

my $LCA = abs_path($ENV{LCA_BIN} // './LCA');
plan skip_all => "LCA executable not found: set LCA_BIN" unless defined $LCA && -x $LCA;

my $LINEAGE = 'k__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus %s';
my $HEADER = join("\t", qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen));
my $LEGACY_HEADER = $HEADER =~ s/qlen$/evalue\tbitscore/r;

sub hit {
	my %h = (query => 'q', subject => 'a', identity => 99, length => 100, qlen => 100, @_);
	my $qend = $h{qend} // $h{length};
	return join("\t", @h{qw(query subject identity length)}, 0, 0, 1, $qend, 1, $h{length}, $h{qlen}) . "\n";
}

sub spew {
	my ($path, $content) = @_;
	open my $fh, '>:raw', $path or die "$path: $!";
	print {$fh} $content;
	close $fh or die "$path: $!";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<:raw', $path or die "$path: $!";
	local $/;
	return scalar <$fh>;
}

# A fresh directory with one hit of query q to subject a, and taxonomy for subjects a and b.
sub fixture {
	my $dir = tempdir('lca-regression-XXXXXX', TMPDIR => 1, CLEANUP => 1);
	my $f = {dir => $dir, mapping => "$dir/input.m8", taxonomy => "$dir/reference.tax", output => "$dir/output.tsv"};
	spew($f->{mapping}, hit());
	spew($f->{taxonomy}, "a\t" . sprintf($LINEAGE, 'alpha') . "\nb\t" . sprintf($LINEAGE, 'beta') . "\n");
	return $f;
}

# Runs LCA and checks its exit status; returns stdout and stderr.
sub run_lca {
	my ($f, %o) = @_;
	my $expected = $o{expected} // 0;
	my @cmd = ($LCA, '-i', $o{mapping} // $f->{mapping}, '-r', $o{taxonomy} // $f->{taxonomy},
		'-o', $o{output} // $f->{output}, @{$o{args} // []});
	my $pid = fork() // die "fork: $!";
	if ($pid == 0) {
		open STDOUT, '>', "$f->{dir}/stdout" or die $!;
		open STDERR, '>', "$f->{dir}/stderr" or die $!;
		exec {$cmd[0]} @cmd or exit 127;
	}
	waitpid($pid, 0);
	my $status = $? >> 8;
	my ($out, $err) = (slurp("$f->{dir}/stdout"), slurp("$f->{dir}/stderr"));
	is($status, $expected, "exit status $expected: @cmd[7 .. $#cmd]") or diag($out . $err);
	return ($out, $err);
}

sub rows {
	my ($f) = @_;
	return [map { [split /\t/, $_, -1] } split /\n/, slurp($f->{output})];
}

sub tail_fields {
	my ($row, $count) = @_;
	return [@{$row}[-$count .. -1]];
}

subtest 'README quick example' => sub {
	my $f = fixture();
	spew($f->{taxonomy}, "refA\t" . sprintf($LINEAGE, 'alpha') . "\nrefB\t" . sprintf($LINEAGE, 'beta') . "\n");
	spew($f->{mapping}, "$HEADER\n" . hit(query => 'query1', subject => 'refA')
		. hit(query => 'query2', subject => 'refA') . hit(query => 'query2', subject => 'refB')
		. "query3\trefA\t96\t100\t4\t0\t1\t100\t1\t100\t100\n");
	run_lca($f, args => ['-showHitRead', '-reportID']);
	is(slurp($f->{output}), join('', map { join("\t", @$_) . "\n" }
		[qw(OTU Domain Phylum Class Order Family Genus Species Hit2DB %ID)],
		['query1', qw(Bacteria P C O F Genus), 'Genus alpha', 'refA', '99.000000'],
		['query2', qw(Bacteria P C O F Genus ? ? 99.000000)],
		['query3', qw(Bacteria P C O F Genus ? refA 96.000000)]), 'documented output');
};

subtest 'tab header, numeric padding and CRLF' => sub {
	my $f = fixture();
	my @fields = split /\t/, hit() =~ s/\n$//r;
	$_ = " $_ " for @fields[2 .. $#fields];
	spew($f->{mapping}, "  # comment\r\n\t \r\n$HEADER\r\n" . join("\t", @fields) . "\r\n");
	run_lca($f);
	is_deeply(rows($f)->[1], ['q', qw(Bacteria P C O F Genus), 'Genus alpha'], 'padded CRLF row parsed');
};

subtest 'legacy layout skips coverage filtering' => sub {
	my $f = fixture();
	spew($f->{mapping}, "$LEGACY_HEADER\nq\ta\t99\t100\t0\t0\t101\t200\t1\t100\t1e-35\t152\n");
	run_lca($f, args => ['-cover', 1]);
	is(rows($f)->[1][-1], 'Genus alpha', 'legacy row assigned');
};

subtest 'malformed rows are skipped with their line number' => sub {
	my $f = fixture();
	for my $invalid (hit() =~ s/99/nan/r, "q a 99\n", "q\ta\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t80\n",
		"q\ta\t99\t100\t0\t0\t0\t100\t1\t100\t100\n", "nohit\t*\t0.0\t0\t0\t0\t0\t0\t0\t0\t100\n") {
		spew($f->{mapping}, hit(query => 'first') . $invalid . hit(query => 'last'));
		my (undef, $err) = run_lca($f);
		like($err, qr/line 2\b/, 'warning names the skipped line');
		like($err, qr/skipped 1 malformed m8 record/, 'summary counts the skipped record');
		is_deeply([map { $_->[0] } @{rows($f)}[1 .. 2]], [qw(first last)], 'valid rows around it are assigned');
	}
};

subtest 'repeated column headers are skipped silently' => sub {
	my $f = fixture();
	spew($f->{mapping}, "$HEADER\n" . hit(query => 'q1') . "$HEADER\n" . hit(query => 'q2'));
	my (undef, $err) = run_lca($f);
	unlike($err, qr/skipping/i, 'no warning for a header');
	is(scalar @{rows($f)}, 3, 'both queries assigned');
};

subtest 'a mapping file without valid records fails' => sub {
	my $f = fixture();
	spew($f->{mapping}, "q a 99 100 0 0 1 100 1 100 100\nq b 99 100 0 0 1 100 1 100 100\n");
	my (undef, $err) = run_lca($f, expected => 25);
	like($err, qr/No valid m8 records.*line 1/, 'error names the first malformed line');
};

subtest 'empty tab-separated fields are not collapsed' => sub {
	my $f = fixture();
	for my $column (0 .. 10) {
		my @fields = split /\t/, hit() =~ s/\n$//r;
		$fields[$column] = '';
		spew($f->{mapping}, join("\t", @fields) . "\n");
		my (undef, $err) = run_lca($f, expected => 25);
		like($err, qr/line 1\b/, "empty column $column reported");
	}
};

subtest 'reported Illumina header preserves query identity' => sub {
	my $f = fixture();
	my $spot = 'LH00409:413:22WTJGLT4:6:2481:9713:19847';
	my $annotated = "$spot 1:N:0:CGTATCTC+CTCGAACA";
	my $other = "$spot 2:N:0:CGTATCTC+CTCGAACA";
	my $subject = 'OBEP010038528';
	spew($f->{taxonomy}, "$subject\t" . sprintf($LINEAGE, 'alpha') . "\n");
	spew($f->{mapping}, hit(query => $spot, subject => $subject, identity => 95.3, length => 150, qlen => 150)
		. "$annotated\t$subject\t 95.3 \t150\t7\t0\t1\t150\t795\t944\t150\n"
		. hit(query => $other, subject => $subject, identity => 95.3, length => 150, qlen => 150));
	run_lca($f, args => ['-showHitRead', '-reportID']);
	my @rows = @{rows($f)}[1 .. 3];
	is_deeply([map { $_->[0] } @rows], [$spot, $annotated, $other], 'query identifiers kept whole');
	is_deeply(tail_fields($_, 2), [$subject, '95.300000'], 'subject and identity') for @rows;
};

subtest 'subject descriptions are preserved' => sub {
	my $f = fixture();
	my ($first, $second) = ('reference description A', 'reference description B');
	spew($f->{taxonomy}, "$first\t" . sprintf($LINEAGE, 'alpha') . "\n$second\t" . sprintf($LINEAGE, 'beta') . "\n");
	spew($f->{mapping}, hit(subject => $first) . hit(subject => $second));
	run_lca($f, args => ['-showHitRead']);
	is_deeply(tail_fields(rows($f)->[1], 2), ['?', '?'], 'consensus of two subjects');
	run_lca($f, args => ['-reportBestHit', '-showHitRead']);
	is_deeply(tail_fields(rows($f)->[1], 2), ['?', $first], 'tied best hits fall back to genus');
};

subtest 'legacy identifiers with spaces' => sub {
	my $f = fixture();
	my ($query, $subject) = ('query with barcode', 'reference with description');
	spew($f->{taxonomy}, "$subject\t" . sprintf($LINEAGE, 'alpha') . "\n");
	spew($f->{mapping}, "$LEGACY_HEADER\n$query\t$subject\t99\t100\t0\t0\t1\t100\t1\t100\t 1e-35 \t 152 \n");
	run_lca($f, args => ['-showHitRead']);
	is(rows($f)->[1][0], $query, 'query kept');
	is_deeply(tail_fields(rows($f)->[1], 2), ['Genus alpha', $subject], 'subject kept');
};

subtest 'non-contiguous queries fail' => sub {
	my $f = fixture();
	spew($f->{mapping}, hit(query => 'q1') . hit(query => 'q2') . hit(query => 'q1'));
	run_lca($f, expected => 27);
};

subtest 'non-contiguous detection covers the last 100,000 completed queries' => sub {
	my $f = fixture();
	my $others = sub { join('', map { hit(query => "other$_") } 1 .. $_[0]) };
	spew($f->{mapping}, hit(query => 'q1') . $others->(99_999) . hit(query => 'q1'));
	run_lca($f, expected => 27);
	spew($f->{mapping}, hit(query => 'q1') . $others->(100_000) . hit(query => 'q1'));
	run_lca($f);
	is(scalar(grep { $_->[0] eq 'q1' } @{rows($f)}), 2, 'repeat beyond the window is not detected');
};

subtest 'tied best hits combine by consensus' => sub {
	my $f = fixture();
	for my $order (['a', 'b'], ['b', 'a']) {
		for my $filtering ([], ['-no_bl_filter']) {
			spew($f->{mapping}, join('', map { hit(subject => $_) } @$order));
			run_lca($f, args => ['-reportBestHit', '-showHitRead', @$filtering]);
			is_deeply(tail_fields(rows($f)->[1], 2), ['?', 'a'], "order @$order @$filtering");
		}
	}
	spew($f->{taxonomy}, "a0\t" . sprintf($LINEAGE, 'beta') . "\n"
		. join('', map { "b$_\t" . sprintf($LINEAGE, 'alpha') . "\n" } 1 .. 4));
	spew($f->{mapping}, join('', map { hit(subject => $_) } qw(b3 a0 b1 b4 b2)));
	run_lca($f, args => ['-reportBestHit', '-showHitRead']);
	is_deeply(tail_fields(rows($f)->[1], 2), ['?', 'a0'], 'disagreeing species stop at genus');
	run_lca($f, args => ['-reportBestHit', '-showHitRead', '-LCAfrac', '0.8']);
	is_deeply(tail_fields(rows($f)->[1], 2), ['Genus alpha', 'b1'], 'Hit2DB names a subject of the consensus');
	spew($f->{mapping}, hit(subject => 'a0', qend => 90) . hit(subject => 'b2'));
	run_lca($f, args => ['-reportBestHit', '-showHitRead']);
	is_deeply(tail_fields(rows($f)->[1], 2), ['Genus alpha', 'b2'], 'higher query coverage wins at equal identity');
};

subtest 'ineligible duplicates do not hide valid hits' => sub {
	my $f = fixture();
	for my $ineligible (hit(identity => 100, length => 50), hit(identity => 100, qlen => 300)) {
		for my $content ($ineligible . hit(), hit() . $ineligible) {
			spew($f->{mapping}, $content);
			run_lca($f, args => ['-reportID']);
			is_deeply(tail_fields(rows($f)->[1], 2), ['Genus alpha', '99.000000'], 'eligible duplicate used');
		}
	}
};

sub permutations {
	my @items = @_;
	return [@items] if @items <= 1;
	return map { my $i = $_; map { [$items[$i], @$_] } permutations(@items[grep { $_ != $i } 0 .. $#items]) } 0 .. $#items;
}

subtest 'relative coverage filter precedes deduplication' => sub {
	my $f = fixture();
	my @alignments = (hit(identity => 97.5), hit(subject => 'b', identity => 97.5, length => 80),
		hit(subject => 'b', identity => 97));
	for my $order (permutations(@alignments)) {
		spew($f->{mapping}, join('', @$order));
		run_lca($f, args => ['-reportID']);
		is_deeply(tail_fields(rows($f)->[1], 2), ['?', '97.250000'], 'same result in every order');
	}
};

subtest 'highest identity anchors the filter despite lower coverage' => sub {
	my $f = fixture();
	spew($f->{taxonomy}, "close\t" . sprintf($LINEAGE, 'alpha') . "\nfar\tk__Bacteria;p__P;c__C;o__O2;f__F2;g__G2;s__G2 beta\n");
	for my $order ([qw(close far)], [qw(far close)]) {
		my %hits = (close => hit(subject => 'close', identity => 99, length => 200, qlen => 250),
			far => hit(subject => 'far', identity => 90, length => 245, qlen => 250));
		spew($f->{mapping}, join('', @hits{@$order}));
		run_lca($f, args => ['-showHitRead', '-reportID']);
		is_deeply(tail_fields(rows($f)->[1], 3), ['Genus alpha', 'close', '99.000000'], "order @$order");
	}
};

subtest 'duplicate subjects do not gain extra votes' => sub {
	my $f = fixture();
	spew($f->{mapping}, hit() x 20 . hit(subject => 'b'));
	run_lca($f, args => ['-showHitRead']);
	is_deeply(tail_fields(rows($f)->[1], 2), ['?', '?'], 'one vote per subject');
};

subtest 'exact relative coverage threshold is inclusive' => sub {
	my $f = fixture();
	spew($f->{mapping}, hit(identity => 97) . hit(subject => 'b', identity => 97, length => 85));
	run_lca($f);
	is(rows($f)->[1][-1], '?', 'both subjects vote');
};

subtest 'identity cutoffs are inclusive' => sub {
	my $f = fixture();
	spew($f->{mapping}, hit(query => 'exact97', identity => 97) . hit(query => 'exact95', identity => 95));
	run_lca($f);
	is_deeply(tail_fields(rows($f)->[1], 2), ['Genus', 'Genus alpha'], '97.0 reaches species');
	is_deeply(tail_fields(rows($f)->[2], 2), ['Genus', '?'], '95.0 reaches genus');
};

subtest 'multi-database selection uses eligible ranks' => sub {
	my $f = fixture();
	my ($mapping2, $taxonomy2) = ("$f->{dir}/second.m8", "$f->{dir}/second.tax");
	spew($f->{mapping}, hit(identity => 94));
	spew($mapping2, hit(subject => 'c', identity => 96));
	spew($taxonomy2, "c\tk__Bacteria;p__P;c__C;o__O;f__F;g__BetterGenus\n");
	for my $pair ([[$f->{mapping}, $mapping2], [$f->{taxonomy}, $taxonomy2]],
		[[$mapping2, $f->{mapping}], [$taxonomy2, $f->{taxonomy}]]) {
		run_lca($f, args => ['-reportID'], mapping => join(',', @{$pair->[0]}), taxonomy => join(',', @{$pair->[1]}));
		is_deeply(tail_fields(rows($f)->[1], 3), ['BetterGenus', '?', '96.000000'], 'better database selected');
	}
};

subtest 'identity mask matches matrix and hit pattern' => sub {
	my $f = fixture();
	my $pattern = "$f->{dir}/pattern.tsv";
	spew($f->{mapping}, hit(identity => 97) . hit(subject => 'b', identity => 96));
	run_lca($f, args => ['-matHigh', '-reportHitPattern', $pattern]);
	is(rows($f)->[1][-1], '?', 'species masked');
	is(slurp("$f->{output}_Species"), "Bacteria;P;C;O;F;Genus;?\t1\n", 'matrix masked');
	is((split /\t/, (split /\n/, slurp($pattern))[1])[-1], '6', 'hit pattern depth');
};

subtest 'hit pattern identity matches the main output' => sub {
	my $f = fixture();
	my $pattern = "$f->{dir}/pattern.tsv";
	spew($f->{mapping}, hit(identity => 97.123456));
	run_lca($f, args => ['-reportID', '-reportHitPattern', $pattern]);
	is(rows($f)->[1][-1], '97.123456', 'main output');
	is((split /\t/, (split /\n/, slurp($pattern))[1])[1], '97.123456', 'hit pattern');
};

subtest 'decimal identity threshold is inclusive' => sub {
	my $f = fixture();
	spew($f->{taxonomy}, join('', map { "s$_\t" . sprintf($LINEAGE, 'alpha') . "\n" } 0 .. 29));
	for my $count (1, 9, 30) {
		spew($f->{mapping}, join('', map { hit(subject => "s$_", identity => 97.2) } 0 .. $count - 1));
		run_lca($f, args => ['-id', '97.2,95,93,91,88,78,0', '-reportID']);
		is_deeply(tail_fields(rows($f)->[1], 2), ['Genus alpha', '97.200000'], "$count subject(s)");
	}
};

subtest 'missing species is not evidence of a certain species' => sub {
	my $f = fixture();
	spew($f->{taxonomy}, "a\t" . sprintf($LINEAGE, 'sp.') . ";t__strainA\n"
		. "b\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__;t__strainB\n");
	spew($f->{mapping}, hit() . hit(subject => 'b'));
	run_lca($f, args => ['-tdep', 8]);
	is(rows($f)->[1][-2], 'Genus sp.', 'uncertain species kept');
};

subtest 'names below a missing rank are used consistently' => sub {
	my $f = fixture();
	my $gapped = "k__Bacteria;p__P;c__C;f__F;g__G;s__G alpha";
	spew($f->{taxonomy}, "g1\t$gapped\ng2\t$gapped\nh1\tk__Bacteria;p__P;c__C;o__O;f__F2;g__G2;s__G2 beta\n");
	spew($f->{mapping}, hit(query => 'one', subject => 'g1') . hit(query => 'two', subject => 'g1')
		. hit(query => 'two', subject => 'g2') . hit(query => 'conflict', subject => 'g1')
		. hit(query => 'conflict', subject => 'h1'));
	run_lca($f);
	my @rows = @{rows($f)}[1 .. 3];
	my @expected = (qw(Bacteria P C ? F G), 'G alpha');
	is_deeply([@{$rows[0]}[1 .. 7]], \@expected, 'one subject');
	is_deeply([@{$rows[1]}[1 .. 7]], \@expected, 'two subjects with the same lineage');
	is_deeply([@{$rows[2]}[1 .. 7]], [qw(Bacteria P C O ? ? ?)], 'conflicting deeper ranks stop the LCA');
};

subtest 'placeholder species and unknown rank names' => sub {
	my $f = fixture();
	my $family = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae';
	spew($f->{taxonomy}, "x\t$family;g__;s__Lachnospiraceae bacterium 3_1_46FAA\n"
		. "y\t$family;g__Blautia;s__Blautia obeum\nu\t$family;g__unknown;s__Unidentified\n");
	spew($f->{mapping}, hit(query => 'vs_certain', subject => 'x') . hit(query => 'vs_certain', subject => 'y')
		. hit(query => 'placeholder_only', subject => 'x') . hit(query => 'unknown_names', subject => 'u'));
	run_lca($f);
	my @rows = @{rows($f)}[1 .. 3];
	is_deeply(tail_fields($rows[0], 2), ['Blautia', 'Blautia obeum'], 'certain species outvotes a placeholder');
	is_deeply(tail_fields($rows[1], 2), ['?', 'Lachnospiraceae bacterium 3_1_46FAA'], 'placeholder kept untruncated alone');
	is_deeply(tail_fields($rows[2], 2), ['?', '?'], 'unknown and unidentified are missing names');
};

subtest 'taxonomy columns after the lineage are ignored' => sub {
	my $f = fixture();
	spew($f->{taxonomy}, "a\t" . sprintf($LINEAGE, 'alpha') . "\textra_column\nb\t" . sprintf($LINEAGE, 'beta') . "\t\n");
	my (undef, $err) = run_lca($f);
	like($err, qr/ignored columns after the lineage on 1 taxonomy line.*first at line 1/, 'warning');
	is_deeply(rows($f)->[1], ['q', qw(Bacteria P C O F Genus), 'Genus alpha'], 'row keeps eight columns');
};

subtest 'malformed taxonomy lines are skipped' => sub {
	my $f = fixture();
	spew($f->{taxonomy}, "no tab here\n\tmissing identifier\na\t" . sprintf($LINEAGE, 'alpha') . "\n");
	my (undef, $err) = run_lca($f);
	like($err, qr/skipped 2 taxonomy line.*first at line 1/, 'warning');
	is(rows($f)->[1][-1], 'Genus alpha', 'valid entry used');
	spew($f->{taxonomy}, "no tab here\n");
	run_lca($f, expected => 12);
};

subtest 'thread option is accepted silently' => sub {
	my $f = fixture();
	my (undef, $err) = run_lca($f, args => ['-t', 4]);
	is($err, '', 'no message about threads');
	is(rows($f)->[1][-1], 'Genus alpha', 'assignment unchanged');
};

# A taxonomy of 300 references with missing ranks and placeholder species, and 3,000 queries
# with 1-8 hits each; references and hits are split by parity for a two-database run.
sub threaded_fixture {
	my ($f) = @_;
	srand(7);
	my (%tax, %m8);
	for my $r (1 .. 300) {
		my $genus = 'G' . ($r % 12);
		my $species = $r % 9 == 0 ? "$genus sp. X$r" : $r % 5 == 0 ? "$genus bacterium" : "$genus s" . ($r % 4);
		my $order = $r % 17 == 0 ? '' : 'o__O' . ($r % 3) . ';';
		$tax{$r % 2} .= "r$r\tk__Bacteria;p__P" . ($r % 2) . ";c__C;${order}f__F" . ($r % 3) . ";g__$genus;s__$species\n";
	}
	for my $q (1 .. 3000) {
		for (1 .. 1 + int(rand(8))) {
			my $r = 1 + int(rand(300));
			my $length = 90 + int(rand(11));
			$m8{$r % 2}{$q} .= join("\t", "q$q", "r$r", sprintf('%.2f', 95 + rand(5)), $length, 1, 0, 1, $length, 1, $length, 100) . "\n";
		}
	}
	my %paths = (taxonomy => "$f->{dir}/all.tax", mapping => "$f->{dir}/all.m8");
	spew($paths{taxonomy}, $tax{0} . $tax{1});
	spew($paths{mapping}, join('', map { ($m8{0}{$_} // '') . ($m8{1}{$_} // '') } 1 .. 3000));
	for my $parity (0, 1) {
		spew("$f->{dir}/db$parity.tax", $tax{$parity});
		spew("$f->{dir}/db$parity.m8", join('', map { $m8{$parity}{$_} // '' } 1 .. 3000));
	}
	return \%paths;
}

subtest 'thread count does not change any output' => sub {
	my $f = fixture();
	my $paths = threaded_fixture($f);
	my @modes = (['-showHitRead', '-reportID', '-matHigh'], ['-reportBestHit', '-showHitRead', '-reportID'],
		['-no_bl_filter', '-LCAfrac', '0.6', '-reportID']);
	for my $mode (0 .. $#modes) {
		my %outputs;
		for my $threads (1, 3, 4) {
			my ($output, $pattern) = ("$f->{dir}/m${mode}_t$threads.tsv", "$f->{dir}/p${mode}_t$threads.tsv");
			run_lca($f, %$paths, output => $output, args => [@{$modes[$mode]}, '-reportHitPattern', $pattern, '-t', $threads]);
			$outputs{$threads} = join("\n--\n", map { -e $_ ? slurp($_) : '' }
				$output, $pattern, "${output}_Species", "${output}_Hit2DB");
		}
		cmp_ok(length $outputs{1}, '>', 30_000, "mode $mode produced output");
		is($outputs{$_}, $outputs{1}, "mode $mode: -t $_ matches -t 1") for 3, 4;
	}
	my %multi;
	for my $threads (1, 4) {
		my $output = "$f->{dir}/multi_t$threads.tsv";
		run_lca($f, mapping => "$f->{dir}/db0.m8,$f->{dir}/db1.m8", taxonomy => "$f->{dir}/db0.tax,$f->{dir}/db1.tax",
			output => $output, args => ['-showHitRead', '-reportID', '-t', $threads]);
		$multi{$threads} = slurp($output);
	}
	is($multi{4}, $multi{1}, 'multiple databases: -t 4 matches -t 1');
};

subtest 'fatal errors are reported with threads' => sub {
	my $f = fixture();
	my $queries = join('', map { hit(query => "q$_") } 1 .. 3000);
	spew($f->{mapping}, $queries =~ s/^(q2000\t)a/${1}absent/mr);
	my (undef, $err) = run_lca($f, args => ['-t', 4], expected => 74);
	like($err, qr/Could not find tax for Subject absent/, 'missing subject named');
	spew($f->{mapping}, $queries . hit(query => 'q2999'));
	(undef, $err) = run_lca($f, args => ['-t', 4], expected => 27);
	like($err, qr/'q2999' occurs in multiple non-contiguous blocks/, 'non-contiguous query named');
	spew($f->{mapping}, "q a 99\n" x 10);
	run_lca($f, args => ['-t', 4], expected => 25);
	SKIP: {
		skip 'directory read errors are platform dependent', 2 if $^O eq 'MSWin32';
		run_lca($f, args => ['-t', 4], mapping => $f->{dir}, expected => 28);
		run_lca($f, args => ['-t', 4], taxonomy => $f->{dir}, expected => 13);
	}
};

subtest 'empty and single final-line inputs' => sub {
	my $f = fixture();
	spew($f->{mapping}, '');
	run_lca($f);
	is(scalar @{rows($f)}, 1, 'header only');
	spew($f->{mapping}, hit() =~ s/\n$//r);
	run_lca($f);
	is(scalar @{rows($f)}, 2, 'final line without newline');
};

SKIP: {
	skip 'gzip disabled in this build', 1 if $ENV{LCA_TEST_NO_GZIP};
	subtest 'gzip success and corruption' => sub {
		my $f = fixture();
		my $compressed = "$f->{dir}/input.m8.gz";
		my $data = hit();
		gzip(\$data => \my $good) or die $GzipError;
		spew($compressed, $good);
		run_lca($f, mapping => $compressed);
		is(scalar @{rows($f)}, 2, 'compressed input read');
		for my $broken (substr($good, 0, -8), substr($good, 0, -8) . "\0" x 8) {
			spew($compressed, $broken);
			run_lca($f, mapping => $compressed, expected => 28);
		}
	};
}

SKIP: {
	skip 'directory read errors are platform dependent', 1 if $^O eq 'MSWin32';
	subtest 'read errors fail' => sub {
		my $f = fixture();
		run_lca($f, mapping => $f->{dir}, expected => 28);
		run_lca($f, taxonomy => $f->{dir}, expected => 13);
	};
}

subtest 'missing subject fails' => sub {
	my $f = fixture();
	spew($f->{mapping}, hit(subject => 'absent'));
	run_lca($f, expected => 74);
};

subtest 'overlapping outputs preserve files' => sub {
	my $f = fixture();
	my $matrix = "$f->{output}_Domain";
	spew($matrix, hit());
	spew($f->{output}, "existing output\n");
	my %original = map { $_ => slurp($_) } $f->{mapping}, $f->{taxonomy}, $f->{output}, $matrix;
	for my $case ([[], output => $f->{mapping}], [[], output => $f->{taxonomy}],
		[['-reportHitPattern', $f->{mapping}]], [['-reportHitPattern', $f->{output}]],
		[['-matHigh'], mapping => $matrix], [['-matHigh', '-reportHitPattern', $matrix]]) {
		my ($args, %paths) = @$case;
		run_lca($f, args => $args, %paths, expected => 29);
		is(slurp($_), $original{$_}, "$_ unchanged") for sort keys %original;
	}
};

SKIP: {
	skip 'link creation needs platform support', 1 if $^O eq 'MSWin32';
	subtest 'input aliases are protected' => sub {
		my $f = fixture();
		my $original = slurp($f->{mapping});
		symlink($f->{mapping}, "$f->{dir}/alias.m8") or die $!;
		link($f->{mapping}, "$f->{dir}/hardlink.m8") or die $!;
		for my $alias ("$f->{dir}/alias.m8", "$f->{dir}/hardlink.m8", "$f->{dir}/./input.m8") {
			run_lca($f, output => $alias, expected => 29);
			is(slurp($f->{mapping}), $original, "$alias did not overwrite input");
		}
	};
}

SKIP: {
	skip 'requires /dev/full', 1 unless -e '/dev/full';
	subtest 'write errors fail' => sub {
		my $f = fixture();
		run_lca($f, output => '/dev/full', expected => 32);
		run_lca($f, args => ['-reportHitPattern', '/dev/full'], expected => 33);
		symlink('/dev/full', "$f->{output}_Domain") or die $!;
		run_lca($f, args => ['-matHigh'], expected => 34);
	};
}

done_testing();
