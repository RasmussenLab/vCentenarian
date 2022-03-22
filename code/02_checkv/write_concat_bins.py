
#!/bin/python
import os
import vamb
from optparse import OptionParser

### From PHAMB repository 

def write_concat_bins(directory, bins, fastadict, compressed=False, maxbins=250, minsize=5000): 
    """Writes bins as FASTA files in a directory, one file per bin.
    Inputs:
        directory: Directory to create or put files in
        bins: {'name': {set of contignames}} dictionary (can be loaded from
        clusters.tsv using vamb.cluster.read_clusters)
        fastadict: {contigname: FastaEntry} dict as made by `loadfasta`
        compressed: Sequences in dict are compressed [False]
        maxbins: None or else raise an error if trying to make more bins than this [250]
        minsize: Minimum number of nucleotides in cluster to be output [0]
    Output: None
    """
    import os as _os
    import gzip as _gzip
    from vamb.vambtools import FastaEntry
    import random

    # Safety measure so someone doesn't accidentally make 50000 tiny bins
    # If you do this on a compute cluster it can grind the entire cluster to
    # a halt and piss people off like you wouldn't believe.
    if maxbins is not None and len(bins) > maxbins:
        raise ValueError('{} bins exceed maxbins of {}'.format(len(bins), maxbins))

    # Check that the directory is not a non-directory file,
    # and that its parent directory indeed exists
    abspath = _os.path.abspath(directory)
    parentdir = _os.path.dirname(abspath)

    if parentdir != '' and not _os.path.isdir(parentdir):
        raise NotADirectoryError(parentdir)

    if _os.path.isfile(abspath):
        raise NotADirectoryError(abspath)

    if minsize < 0:
        raise ValueError("Minsize must be nonnegative")

    # Check that all contigs in all bins are in the fastadict
    allcontigs = set()

    for contigs in bins.values():
        allcontigs.update(set(contigs))

    allcontigs -= fastadict.keys()
    if allcontigs:
        nmissing = len(allcontigs)
        raise IndexError('{} contigs in bins missing from fastadict'.format(nmissing))

    # Make the directory if it does not exist - if it does, do nothing
    try:
        _os.mkdir(directory)
    except FileExistsError:
        pass
    except:
        raise
    
    bins_entries = []
    # Now actually print all the contigs to files
    for binname, contigs in bins.items():
        
        # Concatenate sequences of the bin
        concat_sequence = bytearray()
        for contig in contigs:
            entry = fastadict[contig]
            if compressed:
                uncompressed = bytearray(_gzip.decompress(entry.sequence))
                concat_sequence += uncompressed
            else:
                uncompressed = bytearray(entry.sequence)
                concat_sequence += uncompressed

        bin_entry = FastaEntry(binname, concat_sequence)           
        # Skip bin if it's too small
        if len(bin_entry.sequence) < minsize:
            continue
        bins_entries.append(bin_entry)

    # Print bin to file
    random.shuffle(bins_entries)
    print('Writing:',len(bins_entries) ,'bins to file')
    filename = _os.path.join(directory, 'darkmatter_bins.1.fna')
    i = 1
    j = 1
    file = open(filename,'w')
    for entry in bins_entries:
        if i % 100000 == 0:
            j += 1 
            file.close()
            filename = _os.path.join(directory, 'darkmatter_bins.' + str(j) + '.fna')
            file = open(filename,'w')
        i += 1
        print(entry.format(),file=file)


def filterclusters(clusters, lengthof,minsize = 2500):
    filtered_bins = dict()
    for medoid, contigs in clusters.items():
        binsize = sum(lengthof[contig] for contig in contigs)

        ### Dark-matter specific
        if binsize < 1000000 and binsize >= minsize:
            filtered_bins[medoid] = contigs
    
    return filtered_bins







def main():

  parser = OptionParser()
  parser.add_option("-f", "--input.contigs", dest="in_contigs", help="Combined contigs fasta file", metavar="FILENAME")
  parser.add_option("-o", "--out.dir", dest="out_directory", help="Directory for output files", metavar="FILENAME")
  parser.add_option("-c", "--input.clusters", dest="in_clusterfile", help="VAMB cluster-file", metavar="FILENAME")
  parser.add_option("-s", "--split", dest="split_character", help="Character for splitting bins", metavar="STRING",default=None)
  (options, args) = parser.parse_args()

  contigs= options.in_contigs
  clusterfile = options.in_clusterfile
  bindir = options.out_directory
  splitcharacter = options.split_character
 
    
  ### Read in contigs clustered by VAMB
  keptcontigs = set()
  with open(clusterfile,'r') as infile:
        for line in infile:
            keptcontigs.add( line.strip().split('\t')[1] )

  ### 
  with vamb.vambtools.Reader(contigs, 'rb') as infile:
        fastadict = vamb.vambtools.loadfasta(infile, keep=keptcontigs)
  
  with vamb.vambtools.Reader(contigs, 'rb') as contigfile:
        tnfs, contignames, contiglengths = vamb.parsecontigs.read_contigs(contigfile)
  lengthof = dict(zip(contignames, contiglengths))

  with open(clusterfile,'r') as infile:
        clusters = vamb.vambtools.read_clusters(infile)
  
  if splitcharacter is None:
    darkmatter_bins = filterclusters(clusters, lengthof)
  else:
    darkmatter_bins =  filterclusters(vamb.vambtools.binsplit(clusters, splitcharacter), lengthof)
  #darkmatter_bins = filterclusters(clusters, lengthof)

  ### Write bins to directory out
  write_concat_bins(bindir, darkmatter_bins, fastadict, compressed=False, maxbins=len(darkmatter_bins), minsize=2500)

if __name__ == '__main__':
  main()