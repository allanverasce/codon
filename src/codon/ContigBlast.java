package codon;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import codon.blast.Blaster;
import codon.blast.BlastingManager;
import codon.blast.BlastingManagerListener;
import codon.blast.BlastingThread;
import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.recognizer.OrfRecognizer;

public class ContigBlast {
	
	ArrayList <String> contigs = new ArrayList <> ();
	BlastingManager manager;
	OrfRecognizer reco;
	MemoryMXBean memoryMXBean = ManagementFactory.getMemoryMXBean ();
	
	public ContigBlast (String file) {
		String workspace = CodonConfiguration.workspace;
		try {
			BufferedReader reader = new BufferedReader (new FileReader (file));
			String line = null;
			while ((line = reader.readLine ()) != null) {
				if (line.charAt (0) != '>') contigs.add (line);
			}
			reader.close ();
			contigs.sort (new Comparator <String> () {
				public int compare (String o1, String o2) {
					if (o1.length () < o2.length ()) return -1;
					if (o1.length () > o2.length ()) return 1;
					return 0;
				}
			});
			File f = new File (file);
			String dirB = f.getName ().substring (0, f.getName ().indexOf ('.'));
			CodonConfiguration.workspace = workspace + "/" + dirB + "/";
			processNextContig ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	int totalAccurateORF = 0;
	int totalBlasted = 0;
	int totalIgnored = 0;
		
	public void processNextContig () {
		if (contigs.size () == 0) return;
		int i =contigs.size () - 1;
		String contig = contigs.remove (i);
		
		OrfRecognizer reco = new OrfRecognizer ();
		reco.loadFastaString (i + "", contig);
		CodonConfiguration.out.println ("======================================================================================");
		CodonConfiguration.out.println ("Processing contig " + i + " from " + contigs.size () + ", contig size: " + contig.length () + ", orfs: " + reco.getOrfCount ());
		CodonConfiguration.out.println ("======================================================================================");
		manager = new BlastingManager (reco, 10, false);
		manager.addBlastingManagerListener (new BlastingManagerListener () {
			public void threadStatustChanged (BlastingManager manager, BlastingThread endedThread) {
				CodonConfiguration.out.println ("Connection status changed: " + endedThread.getId () + " connection = " + endedThread.withoutConnection);
			}
			public void threadEnded (BlastingManager manager, BlastingThread endedThread) {
				ORF orf = endedThread.orf;
				CodonConfiguration.out.println ("-------------------------------------------------------------------------------------------------");
				CodonConfiguration.out.println (Thread.activeCount () + " running threads - " + Blaster.ongoingBlast + " ongoing blasts - " + Blaster.ongoingEntry + " entry requests - Memory usage:" + String.format (": %.2f GB", (double) memoryMXBean.getHeapMemoryUsage ().getUsed () / 1073741824));
				if  (orf.getAccuracy () > 80) CodonConfiguration.out.print ("* ");
				CodonConfiguration.out.println ("BLASTED: " + orf.id + ", Acc:" + orf.getAccuracy () + ", "  + (orf.product.isGene () ? orf.product.gene + "|" : "") + orf.product.productName);
				CodonConfiguration.out.println ("Done: " + manager.done + ", ignored: " + manager.ignored + ", stopped: " + manager.stopped + ", Previously Done: " + manager.previouslyDone + ", to process: " + (manager.toProcess () + manager.ongoing ()));
			}
			public void blastingFinalized (BlastingManager manager) {
				totalBlasted += manager.previouslyDone + manager.done;
				totalIgnored += manager.ignored;
				reco.removeLowAccuracy (80);
				totalAccurateORF += reco.getOrfCount ();
				CodonConfiguration.out.println ("======================================================================================");
				CodonConfiguration.out.println (">>>> Accurate count: " + reco.getOrfCount () + " / " + totalAccurateORF);
				CodonConfiguration.out.println (">>>> Total blasted: " + totalBlasted);
				CodonConfiguration.out.println (">>>> Total ignored: " + totalIgnored);
				processNextContig ();
			}
			public void threadStarted (BlastingManager manager, BlastingThread endedThread) {}
			public void threadRestarted (BlastingManager manager, BlastingThread endedThread) {}
			public void orfIgnored (BlastingManager manager, ORF orf) {}
			public void orfPreviouslyDone (BlastingManager manager, ORF orf, int toProcess) {}
			public void orfRebounded (BlastingManager manager, ORF orf, int start, int stop) {}
			public void orfReboundedAndRename (BlastingManager manager, ORF orf, int start, int stop) {}
			public void orfReincluded(BlastingManager manager, ORF orf) {}
		});
		ArrayList <ORF> orfsToBlast = new ArrayList <ORF> ();
		for (int frame = 0; frame < 6; frame++) {
			List <ORF> orfs = reco.orfsPorFrame [frame];
			for (ORF orf: orfs) {
				if (orf.isRemoved ()) continue;
				if (orf.length < CodonConfiguration.min_orf_length) {
					System.err.println ("Why are you here???");
					continue;
				}
				orfsToBlast.add (orf);
			}
		}
		manager.start (orfsToBlast);
	}

	public static void main (String [] args) {
		CodonConfiguration.out = System.out;
		System.setOut (new PrintStream (new OutputStream () {
			public void write (int arg0) throws IOException {}
		}));
		new ContigBlast ("workspace/META_GENOMA.contig");
	}
	
}
