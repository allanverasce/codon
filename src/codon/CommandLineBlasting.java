package codon;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.util.ArrayList;
import java.util.List;

import codon.blast.Blaster;
import codon.blast.BlastingManager;
import codon.blast.BlastingManagerListener;
import codon.blast.BlastingThread;
import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.recognizer.OrfRecognizer;

public class CommandLineBlasting {
	
	public CommandLineBlasting (String enrtyFile, int threadNumber, boolean upgradeWithSwissProt, boolean fullblast, boolean retreiveRNA, int start, int end) {
		this ("./workspace", enrtyFile, threadNumber, upgradeWithSwissProt, fullblast, retreiveRNA, start, end);		
	}
	
	public CommandLineBlasting (String enrtyFile, int threadNumber) {
		this ("./workspace", enrtyFile, threadNumber, false, false, true, -1, -1);
	}
	
	BlastingManager manager;
	OrfRecognizer reco;
	int start = -1;
	int end = -1;
	boolean fullblast = false;
	boolean upgradeWithSwissProt = false;
	boolean retreiveRNA = true;
	
	boolean lastAcceleratedWithSwissProt = false;
	
	static String exportEmblFile = null;
	static String exportCodonFile = null;
	static int accuracyFilter = -1;
	static int overlapThreeshold = -1;
	static boolean optimizeIntegenicRegion = false;
	
	MemoryMXBean memoryMXBean = ManagementFactory.getMemoryMXBean ();
	
	public CommandLineBlasting (String workspace, String enrtyFile, int threadNumber, boolean accelerateWithSwissProt, boolean fullblast, boolean retreiveRNA, int start, int end) {
		Blaster.acceleratedWithSwissProt = accelerateWithSwissProt;
		lastAcceleratedWithSwissProt = accelerateWithSwissProt;
		CodonConfiguration.workspace = workspace + "/";
		new File (workspace).mkdirs ();
		reco = new OrfRecognizer ();
		if (enrtyFile.contains (".fa")) reco.load (enrtyFile);
		if (enrtyFile.contains (".embl")) reco.loadFromEmbl (new File (enrtyFile));
		if (enrtyFile.contains (".gb")) reco.loadFromEmbl (new File (enrtyFile));
		reco.loadPreannotated ();
		manager = new BlastingManager (reco, threadNumber, fullblast);
		manager.addBlastingManagerListener (new BlastingManagerListener () {
			public void threadStatustChanged (BlastingManager manager, BlastingThread endedThread) {
				CodonConfiguration.out.println ("<STATUS>");
				CodonConfiguration.out.println (endedThread.orf.id + ";" + endedThread.withoutConnection);
				CodonConfiguration.out.flush ();
			}
			public void threadEnded (BlastingManager manager, BlastingThread endedThread) {
				ORF orf = endedThread.orf;
				CodonConfiguration.out.println ("<ENDED>");
				CodonConfiguration.out.println (System.currentTimeMillis () + ";" + orf.id + ";Acc=" + orf.getAccuracy () + ";"  + (orf.product.isGene () ? orf.product.gene + "|" : "") + orf.product.productName);
				CodonConfiguration.out.println ("Done=" + manager.done + ";Ignored=" + manager.ignored + ";Stopped=" + manager.stopped + ";Previously-Done=" + manager.previouslyDone + ";To-Process=" + manager.toProcess ());
				CodonConfiguration.out.println ("Running-threads=" + Thread.activeCount () + ";Ongoing-Blasts=" + Blaster.ongoingBlast + ";Ongoing-Entries=" + Blaster.ongoingEntry
						+ ";Memory-Load=" + String.format ("%.2fGB", (double) memoryMXBean.getHeapMemoryUsage ().getUsed () / 1073741824)
						+ ";Memory-Comitted=" + ((long) (100 * (memoryMXBean.getHeapMemoryUsage ().getUsed () / (double) memoryMXBean.getHeapMemoryUsage ().getCommitted ()))
						+ "%;Memory-Max=" + ((long) (100 * (memoryMXBean.getHeapMemoryUsage ().getUsed () / (double) memoryMXBean.getHeapMemoryUsage ().getMax ())))) + "%");
				if (lastAcceleratedWithSwissProt != Blaster.acceleratedWithSwissProt) {
					CodonConfiguration.out.println ("<SWISS>");
					CodonConfiguration.out.println (Blaster.acceleratedWithSwissProt);
					lastAcceleratedWithSwissProt = Blaster.acceleratedWithSwissProt;
				}
			}
			public void blastingFinalized (BlastingManager manager) {
				CodonConfiguration.out.println ("<FINALIZING>");
				if (exportEmblFile == null && exportCodonFile == null) {
					//System.exit (0);
				}
				else {
					if (accuracyFilter > 0 && accuracyFilter < 100) {
		//					CodonConfiguration.out.println ("Removing orf with accuracy lower than: " + accuracyFilter);
						reco.removeLowAccuracy (accuracyFilter);
					}
					if (overlapThreeshold >= 0 && overlapThreeshold < 100) {
		//					CodonConfiguration.out.println ("Removing overlaps (overlap %): " + overlapThreeshold);
						reco.removeOverlaps (overlapThreeshold);
					}
					if (optimizeIntegenicRegion) {
		//					CodonConfiguration.out.println ("Minimizing intergenic regions");
						reco.selectBetterEntriesFulfillIntergenic (2);
					}
					if (exportEmblFile != null) {
		//					CodonConfiguration.out.println ("Exporting as EMBL " + exportEmblFile);
						reco.exportAsEmbl (new File (exportEmblFile));
					}
					if (exportCodonFile != null) {
		//					CodonConfiguration.out.println ("Exporting as CODON " + exportCodonFile);
						reco.saveAsCodonProject (new File (exportCodonFile));
					}
				}
				CodonConfiguration.out.println ("<FINISHED>");
//				try {
//					Thread.sleep (5000);
//				} catch (InterruptedException e) {
//					e.printStackTrace ();
//				}
//				System.exit (0);
			}
			public void threadStarted (BlastingManager manager, BlastingThread endedThread) {
				CodonConfiguration.out.println ("<STARTED>");
				CodonConfiguration.out.println (endedThread.orf.id);
			}
			public void threadRestarted (BlastingManager manager, BlastingThread endedThread) {
				CodonConfiguration.out.println ("<RESTARTED>");
				CodonConfiguration.out.println (endedThread.orf.id);
			}
			public void orfIgnored (BlastingManager manager, ORF orf) {
				CodonConfiguration.out.println ("<I>");
				CodonConfiguration.out.println (orf.id);
			}
			public void orfPreviouslyDone (BlastingManager manager, ORF orf, int toProcess) {
				CodonConfiguration.out.println ("<P>");
				CodonConfiguration.out.println (orf.id + ";" + toProcess);
			}
			public void orfRebounded(BlastingManager manager, ORF orf, int start, int stop) {
				CodonConfiguration.out.println ("<R>");
				CodonConfiguration.out.println (orf.id + ";" + start + ";" + stop);
			}
			public void orfReboundedAndRename (BlastingManager manager, ORF orf, int start, int stop) {
				CodonConfiguration.out.println ("<N>");
				CodonConfiguration.out.println (orf.id + ";" + start + ";" + stop);
			}
			public void orfReincluded (BlastingManager manager, ORF orf) {
				CodonConfiguration.out.println ("<RI>");
				CodonConfiguration.out.println (orf.id);
			}
		});
	}
	
	public void start () {
		if (retreiveRNA) {
			reco.annotateRNAOrf (false);
//			reco.addTRNAOrf (false);
		}
		ArrayList <ORF> orfsToBlast = new ArrayList <ORF> ();
		for (int frame = 0; frame < 6; frame++) {
			List <ORF> orfs = reco.orfsPorFrame [frame];
			for (ORF orf: orfs) {
				if (orf.isRemoved ()) continue;
				if (end > start && orf.start < start && orf.start >= end) continue;
				orfsToBlast.add (orf);
			}
		}
		manager.start (orfsToBlast);
		BufferedReader reader = new BufferedReader (new InputStreamReader (System.in));
		String line = null;
		try {
			while ((line = reader.readLine ()) != null) {
				if (line.equals ("DIE")) System.exit (0);
				else if (line.indexOf ("TRIES=") == 0) {
					int nb = Integer.parseInt (line.substring (6));
					nb = Math.min (nb, 20);
					manager.setMaxThread (nb);
				}
				else if (line.indexOf ("SWISS=") == 0) {
					Blaster.acceleratedWithSwissProt = Boolean.parseBoolean (line.substring (6));
				}
				else if (line.indexOf ("RESTART=") == 0) {
					String orfId = line.substring (8);
					manager.restart (orfId);
					//Blaster.acceleratedWithSwissProt = Boolean.parseBoolean (line.substring (6));
				}
			}
		} catch (IOException e) {
			e.printStackTrace ();
		}
		System.exit (0);
	}
	
	public static void main (String [] args) {
		CodonConfiguration.ui = false;
		String workspace = "./workspace";
		int threadNumber = 1;
		boolean fullBlast = false;
		boolean accelerateWithSwiss = false;
		boolean retreiveRNA = true;
		int start = -1;
		int end = -1;
		
		for (int i = 0; i < args.length; i++) {
			if (args [i].equalsIgnoreCase ("--fullblast")) fullBlast = true;
			if (args [i].equalsIgnoreCase ("--ignore-rna")) retreiveRNA = false;
			if (args [i].equalsIgnoreCase ("--accelerate-swiss")) accelerateWithSwiss = true;
			if (args [i].indexOf ("-simultanous-blast=") == 0) threadNumber = safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-start=") == 0) start = safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-end=") == 0) end = safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-workspace=") == 0) workspace = args [i].substring (args[i].indexOf ("=") + 1);
			if (args [i].indexOf ("-export-emb=") == 0) exportEmblFile = args [i].substring (args[i].indexOf ("=") + 1).replaceAll ("%20", " ");
			if (args [i].indexOf ("-export-codon=") == 0) exportCodonFile = args [i].substring (args[i].indexOf ("=") + 1).replaceAll ("%20", " ");
			if (args [i].indexOf ("-overlaps-threehold=") == 0) overlapThreeshold = safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-low-accuracy-threehold=") == 0) accuracyFilter = safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-accuracy-tolerance=") == 0) CodonConfiguration.tolerance = safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-accuracy-removing-overlap-during-blast=") == 0) CodonConfiguration.accuracy_level_to_remove_overlap_during_blast = safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("--minimize-intergenic") == 0) optimizeIntegenicRegion = true;
		}
		
		if (args.length == 0) {
			System.err.println ("The command line must end by an existing fasta file or embl file");
			System.exit (0);
		}
		String file = args [args.length - 1];
		file = file.replaceAll ("%20", " ");
		workspace = workspace.replaceAll ("%20", " ");
		if (!new File (file).exists ()) {
			System.err.println ("The file " + file + " does not exist");
			System.exit (0);
		}
		
		CodonConfiguration.out = System.out;
		System.setOut (new PrintStream (new OutputStream () {
			public void write (int arg0) throws IOException {}
		}));
		
		CommandLineBlasting process = new CommandLineBlasting (workspace, file, threadNumber, accelerateWithSwiss, fullBlast, retreiveRNA, start, end);
		process.start ();
	}
	
	static public int safeParamToIntConversion (String v) {
		try {
			return Integer.parseInt (v.substring (v.indexOf ("=") + 1));
		}
		catch (Exception e) {
			System.err.println ("Unvalid " + v + " value");
			System.exit (0);
		}
		return -1;
	}
}
