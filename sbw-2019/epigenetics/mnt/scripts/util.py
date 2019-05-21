#!/usr/bin/env python
# Original https://github.com/JetBrains-Research/washu
# Author oleg.shpynov@jetbrains.com

import getopt
import glob
import os
import re
import subprocess
import sys
import traceback

__author__ = 'oleg.shpynov@jetbrains.com'

help_message = '''
Usage:

python util.py find_input <file>
    Finds input given the file name. Heuristics: among all the files within\
     folder find file with "input" substring and
    most common subsequence with initial file.

python util.py macs_species <genome>
    Converts UCSC genome name to MACS.

python util.py effective_genome_fraction <genome> <chrom.sizes.path>
    Computes effective genome size, required for SICER.
'''


def usage():
    print(help_message)


def lcs(x, y):
    """
    Finds longest common subsequence
    Code adopted from https://en.wikibooks.org/wiki/Algorithm_Implementation/
    Strings/Longest_common_subsequence#Python
    """
    m = len(x)
    n = len(y)
    # An (m+1) times (n+1) matrix
    c = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if x[i - 1] == y[j - 1]:
                c[i][j] = c[i - 1][j - 1] + 1
            else:
                c[i][j] = max(c[i][j - 1], c[i - 1][j])

    def back_track(i, j):
        if i == 0 or j == 0:
            return ""
        elif x[i - 1] == y[j - 1]:
            return back_track(i - 1, j - 1) + x[i - 1]
        else:
            if c[i][j - 1] > c[i - 1][j]:
                return back_track(i, j - 1)
            else:
                return back_track(i - 1, j)

    return len(back_track(m, n))


def is_input(c):
    return re.match('.*input.*', re.sub('.*/', '', str(c)), flags=re.IGNORECASE) is not None


def find_input_name(c, variants):
    if is_input(c):
        return None

    def sort_function(x):
        return lcs(str(c), x.lower())

    inputs = [str(v) for v in variants if is_input(v)]
    if len(inputs) > 0:
        return max(inputs, key=sort_function)
    else:
        return None


def find_input(bam):
    bam_name = os.path.basename(bam).lower()
    if 'input' in bam_name:
        return ''

    # Find all the files within folder
    dirname = os.path.dirname(bam)
    names = [os.path.basename(n) for n in glob.glob('{}/*.bam'.format(dirname))]
    input_name = find_input_name(bam_name, names)
    if input_name is None:
        return ''
    else:
        return input_name


def macs_species(genome):
    """Convert genome to macs2 species encoding"""
    if re.match('^hg[0-9]+$', genome):
        return 'hs'
    elif re.match('^mm[0-9]+$', genome):
        return 'mm'
    raise Exception('Unknown species {}'.format(genome))


def effective_genome_fraction(genome, chrom_sizes_path, pileup_bed):
    """From MACS2 documentation:
    The default hs 2.7e9 is recommended for UCSC human hg18 assembly.
    Here are all precompiled parameters for effective genome size:
    hs: 2.7e9
    mm: 1.87e9
    ce: 9e7
    dm: 1.2e8"""
    chromosomes = re.split('\n',
                           run([['cat', pileup_bed],
                                ['awk', '{print $1}'],
                                ['sort', '--unique']])[0].decode('utf-8').strip())
    chrom_sizes = {}
    with open(chrom_sizes_path, 'r') as f:
        for line in f.readlines():
            chromosome, size = re.split('\t', line.strip())
            chrom_sizes[chromosome] = int(size)

    chromosomes_length = sum([chrom_sizes[c] if c in chrom_sizes else 0 for c in chromosomes])
    genome_length = sum(chrom_sizes.values())

    if genome.startswith('mm'):
        size = 1.87e9
    elif genome.startswith('hg'):
        size = 2.7e9
    else:
        raise Exception('Unknown species {}'.format(genome))
    return (size / genome_length) * (1.0 * chromosomes_length / genome_length)


def run(commands, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """Launches pipe of commands given stdin and final stdout, stderr"""

    # TODO: consider 'plumbum' library instead of this
    # https://plumbum.readthedocs.io/en/latest/local_commands.html

    processes = []
    _stdin = stdin

    for i, cmd in enumerate(commands):
        if i < len(commands) - 1:
            _stdout = subprocess.PIPE
            # Not clear how to collect stderr from chain, let's left
            # None here because result is more consistent:
            # * last cmd stderr is captured
            # * intermediate stderr is not missed, but not captured, goes to
            #    stderr
            # If you feel power, try to fix it + see tests
            _stderr = None
        else:
            _stdout = stdout
            _stderr = stderr

        p = subprocess.Popen(cmd, stdin=_stdin, stdout=_stdout,
                             stderr=_stderr)
        processes.append(p)
        _stdin = p.stdout

    for i in range(0, len(processes)):
        # noinspection PyBroadException
        try:
            if i < len(processes) - 1:
                # Allow p1 to receive a SIGPIPE if p2 exits.
                processes[i].stdout.close()
            else:
                sp = processes[i]
                out, err = sp.communicate()
                return out, err  # , sp.returncode

        except Exception:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print("Error running: {}".format(commands))
            # exc_type below is ignored on 3.5 and later
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                                      # limit=2,
                                      file=sys.stdout)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    # Process help
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            return

    if len(args) == 2 and args[0] == 'find_input':
        print(find_input(args[1]))

    if len(args) == 2 and args[0] == 'macs_species':
        print(macs_species(args[1]))

    if len(args) == 4 and args[0] == 'effective_genome_fraction':
        print(effective_genome_fraction(args[1], args[2], args[3]))


if __name__ == "__main__":
    main()
