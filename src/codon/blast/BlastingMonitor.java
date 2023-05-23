package codon.blast;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Timer;
import java.util.TimerTask;

import codon.CommandLineBlasting;
import codon.data.CodonConfiguration;

public class BlastingMonitor extends Thread {
	
	static boolean standaloneApplication = false;
	
	static String statusFile = "data/blast-monitor.report";
	
	Process proc;
	long lastHit = 0;
	
	int blasting = 0;
	int ended = 0;
	int deadRestart = 0;
	int memoryRestart = 0;
	int timeoutRestart = 0;
	boolean finalized = false;
	String status = "";
	int maxLoad = 0;
	
	int toProcess = 0;
	int ignored = 0;
	int previouslyDone = 0;
	
	boolean firstStart = true;
	
	int toProcessAtStart = 0;
	
	ArrayList <BlastingMonitorListener> listeners = new ArrayList <BlastingMonitorListener> ();
	
	public void addBlastingMonitorListener (BlastingMonitorListener l) {
		listeners.add (l);
	}
	
	public void removeBlastingMonitorListener (BlastingMonitorListener l) {
		listeners.remove (l);
	}
	
	public void fireEnded (String orfId) {
		for (BlastingMonitorListener l: listeners) l.blastEnded (this, orfId); 
	}
	
	public void fireStarted (String orfId) {
		for (BlastingMonitorListener l: listeners) l.blastStarted (this, orfId); 
	}
	
	public void fireRestarted () {
		for (BlastingMonitorListener l: listeners) l.blastRestarted (this); 
	}
	
	public void fireRestarted (String orfId) {
		for (BlastingMonitorListener l: listeners) l.blastRestarted (this, orfId); 
	}
	
	public void fireIgnored (String orfId) {
		for (BlastingMonitorListener l: listeners) l.orfIgnored (this, orfId); 
	}
	
	public void fireReincluded (String orfId) {
		for (BlastingMonitorListener l: listeners) l.orfReincluded (this, orfId); 
	}
	
	public void fireRebounded (String orfId, int start, int stop) {
		for (BlastingMonitorListener l: listeners) l.orfRebounded (this, orfId, start, stop); 
	}
	
	public void fireReboundedAndRenamed (String orfId, int start, int stop) {
		for (BlastingMonitorListener l: listeners) l.orfReboundedAndRenamed (this, orfId, start, stop); 
	}
	
	public void firePreviouslyDone (String orfId) {
		for (BlastingMonitorListener l: listeners) l.blastPreviouslyDone (this, orfId); 
	}
	
	public void fireFinalized () {
		for (BlastingMonitorListener l: listeners) l.blastingFinalized (this); 
	}
	
	
	public void fireStatusChanged (String orfId, boolean status) {
		for (BlastingMonitorListener l: listeners) l.statusChanged (this, orfId, status);
	}

	
	class ErrorControlerReader extends Thread {
		public void run () {
			InputStream inP = proc.getErrorStream ();
			BufferedReader reader = new BufferedReader (new InputStreamReader (inP));
			String line;
			try {
				while ((line = reader.readLine ()) != null) {
					if (blastErrors) System.err.println ("BLAST-ERROR:" + line);
				}
			}
			catch (IOException io) {
			}
		}
	}
	
	public void updateReport () {
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter (statusFile));
			writer.write ("Blasting: " + blasting + " / " + simultaneous);
			writer.newLine ();
			writer.write ("Ended: " + ended);
			writer.newLine ();
			writer.write ("Remaining: " + toProcess);
			writer.newLine ();
			writer.write ("Previously done: " + previouslyDone);
			writer.newLine ();
			writer.write ("Ignored: " + ignored);
			writer.newLine ();
			writer.write ("Max memory load: " + maxLoad);
			writer.newLine ();
			writer.write ("Dead restarts: " + deadRestart);
			writer.newLine ();
			writer.write ("Overload memory restarts: " + memoryRestart);
			writer.newLine ();
			writer.write ("Timeout restarts: " + timeoutRestart);
			writer.newLine ();
			writer.write (status);
			writer.newLine ();
			writer.close ();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	class OutputControlerReader extends Thread {
		public void run () {
			InputStream inP = proc.getInputStream ();
			BufferedReader reader = new BufferedReader (new InputStreamReader (inP));
			String line;
			try {
				while ((line = reader.readLine ()) != null) {
					if (verbose) System.out.println (line);
					if (line.equals ("<STARTED>")) {
						lastHit = System.currentTimeMillis ();
						blasting ++;
						line = reader.readLine (); 
						if (verbose) System.out.println (line);
						fireStarted (line);
						updateReport ();
					}
					else if (line.equals ("<RESTARTED>")) {
						line = reader.readLine (); 
						if (verbose) System.out.println (line);
						fireRestarted (line);
						updateReport ();
					}
					else if (line.equals ("<STATUS>")) {
						line = reader.readLine (); 
						if (verbose) System.out.println (line);
						String s [] = line.split (";");
						fireStatusChanged (s [0], Boolean.parseBoolean (s [1]));
					}
					else if (line.equals ("<P>")) {
						line = reader.readLine (); 
						if (verbose) System.out.println (line);
						String s [] = line.split (";");
						firePreviouslyDone (s [0]);
						previouslyDone ++;
						toProcess = Integer.parseInt (s [1]);
						updateReport ();
					}
					else if (line.equals ("<I>")) {
						line = reader.readLine (); 
						if (verbose) System.out.println (line);
						fireIgnored (line);
						ignored ++;
						toProcess --;
						updateReport ();
					}
					else if (line.equals ("<RI>")) {
						line = reader.readLine (); 
						if (verbose) System.out.println (line);
						fireReincluded (line);
						ignored --;
						toProcess ++;
					}
					else if (line.equals ("<R>")) {
						line = reader.readLine (); 
						if (verbose) System.out.println (line);
						String s [] = line.split (";");
						fireRebounded (s [0], Integer.parseInt (s[1]), Integer.parseInt (s[2]));
					}
					else if (line.equals ("<N>")) {
						line = reader.readLine (); 
						if (verbose) System.out.println (line);
						String s [] = line.split (";");
						fireReboundedAndRenamed (s [0], Integer.parseInt (s[1]), Integer.parseInt (s[2]));
					}
					else if (line.equals ("<ENDED>")) {
						toProcess --;
						lastHit = System.currentTimeMillis ();
						status = "";
						status += reader.readLine () + ";";
						status += reader.readLine () + ";";
						status += reader.readLine ();
						String s [] = status.split (";");
						if (verbose) System.out.println (status);
						for (int i = 0; i < s.length; i++) {
							if (s [i].contains ("Memory-Comitted=")) {
								int load = Integer.parseInt (s [i].substring (16, s [i].indexOf ("%")));
								if (load > maxLoad) maxLoad = load;
								if (load > 95) {
									System.err.println (status);
									System.err.println ("Process memory overloaded: restarting");
									memoryRestart ++;
									kill ();
								}
							}
							if (s [i].contains ("To-Process=")) {
								toProcess = Integer.parseInt (s [i].substring (11));
							}
						}
						blasting --;
						ended ++;
						fireEnded (s [1]);
						updateReport ();
					}
					else if (line.equals ("<FINALIZING>")) {
						finalized = true;
						toProcess = 0;
						timer.cancel ();
						fireFinalized ();
					} 
					else if (line.equals ("<FINISHED>")) {
						kill ();
						if (standaloneApplication) System.exit (0);
					}
					else System.err.println ("OutputControlerReader: CommandLine unexpected response: " + line);
				}
				reader.close ();
			} catch (IOException e) {
				e.printStackTrace ();
			}
		};
	};
	PrintStream  ouP = null;
	
	public void restartOrf (String orfId) {
		ouP.println ("RESTART=" + orfId);
		ouP.flush ();
	}

	public void restartAll () {
		kill ();
	}
	
	public void kill () {
		ouP.println ("DIE");
		ouP.flush ();
		try {
			Thread.sleep (1000);
			if (proc.isAlive ()) {
				Thread.sleep (2000);
				if (proc.isAlive ()) {
					proc.destroy ();
					Thread.sleep (2000);
					if (proc.isAlive ()) {
						proc.destroyForcibly ();
					}
				}
			}
		} catch (InterruptedException e) {
			e.printStackTrace ();
		}
	}
	
	String file = "";
	int simultaneous = 1;
	private boolean accelerateSwiss = true;
	boolean fullBlast = false;
	boolean retreiveRNA = true;
	int start = -1;
	int end = -1;
	
	public void startProcess () {
		maxLoad = 0;
		blasting = 0;
		toProcess = toProcessAtStart;
		ignored = 0;
		previouslyDone = 0;
		ended = 0;
		try {
			String os = System.getProperty ("os.name");
			String wsp = CodonConfiguration.workspace;
			wsp = wsp.replaceAll (" ", "%20");
			file = file.replaceAll (" ", "%20");
			if (os.contains ("Windows")) {
				proc = Runtime.getRuntime ().exec ("java -Xmx2548m "
						+ "-cp \"lib\\activation-1.1.jar;lib\\animal-sniffer-annotations-1.14.jar;lib\\antlr4-4.5.3.jar;lib\\aopalliance-1.0.jar;lib\\aopalliance-repackaged-2.4.0-b10.jar;lib\\avro-1.7.6.jar;lib\\avro-ipc-1.7.6.jar;lib\\axis2-1.6.2.jar;lib\\checker-compat-qual-2.0.0.jar;lib\\classmate-0.8.0.jar;lib\\commons-beanutils-1.9.2.jar;lib\\commons-collections-3.2.1.jar;lib\\commons-compress-1.12.jar;lib\\commons-io-2.4.jar;lib\\commons-lang-2.6.jar;lib\\commons-lang3-3.6.jar;lib\\commons-logging-1.2.jar;lib\\commons-math3-3.6.1.jar;lib\\commons-net-3.6.jar;lib\\commons-text-1.8.jar;lib\\error_prone_annotations-2.1.3.jar;lib\\failsafe-1.0.3.jar;lib\\FastInfoset-1.2.13.jar;lib\\guava-24.1.1-jre.jar;lib\\guice-4.0.jar;lib\\guice-multibindings-4.0.jar;lib\\hibernate-validator-5.0.1.Final.jar;lib\\hk2-api-2.4.0-b10.jar;lib\\hk2-locator-2.4.0-b10.jar;lib\\hk2-utils-2.4.0-b10.jar;lib\\httpclient-4.5.3.jar;lib\\httpcore-4.4.6.jar;lib\\httpmime-4.5.3.jar;lib\\istack-commons-runtime-3.0.5.jar;lib\\j2objc-annotations-1.1.jar;lib\\jackson-annotations-2.9.5.jar;lib\\jackson-core-2.9.5.jar;lib\\jackson-core-asl-1.9.13.jar;lib\\jackson-databind-2.9.4.jar;lib\\jackson-mapper-asl-1.9.13.jar;lib\\japi-1.0.30.jar;lib\\javaparser-1.0.11.jar;lib\\javassist-3.18.1-GA.jar;lib\\javax.annotation-api-1.2.jar;lib\\javax.inject-1.jar;lib\\javax.inject-2.4.0-b10.jar;lib\\javax.mail-1.6.1.jar;lib\\javax.mail-api-1.6.1.jar;lib\\javax.ws.rs-api-2.0.1.jar;lib\\jaxb2-basics-1.11.1.jar;lib\\jaxb2-basics-runtime-1.11.1.jar;lib\\jaxb2-basics-tools-1.11.1.jar;lib\\jaxb-api-2.3.0.jar;lib\\jaxb-core-2.3.0.jar;lib\\jaxb-impl-2.3.0.jar;lib\\jaxb-runtime-2.3.0.jar;lib\\jboss-logging-3.1.1.GA.jar;lib\\jcl-over-slf4j-1.7.7.jar;lib\\jcommander-1.32.jar;lib\\jersey-client-2.17.jar;lib\\jersey-common-2.17.jar;lib\\jersey-guava-2.17.jar;lib\\jersey-media-jaxb-2.17.jar;lib\\jetty-6.1.26.jar;lib\\jetty-util-6.1.26.jar;lib\\jsr305-1.3.9.jar;lib\\log4j-api-2.12.1.jar;lib\\log4j-core-2.12.1.jar;lib\\logback-classic-1.2.3.jar;lib\\logback-core-1.2.3.jar;lib\\metrics-core-3.2.2.jar;lib\\netty-3.4.0.Final.jar;lib\\noggit-0.8.jar;lib\\osgi-resource-locator-1.0.1.jar;lib\\paranamer-2.3.jar;lib\\servlet-api-2.5-20081211.jar;lib\\slf4j-api-1.7.25.jar;lib\\snappy-java-1.0.5.jar;lib\\solr-solrj-7.2.0.jar;lib\\stax2-api-3.1.4.jar;lib\\stax-ex-1.7.8.jar;lib\\txw2-2.3.0.jar;lib\\validation-api-1.1.0.Final.jar;lib\\velocity-1.7.jar;lib\\woodstox-core-asl-4.4.1.jar;lib\\zookeeper-3.4.10.jar;bin\" -Dswing.defaultlaf=javax.swing.plaf.nimbus.NimbusLookAndFeel"
						+ " codon.CommandLineBlasting"
						+ " -simultanous-blast=" + simultaneous
						+ " -workspace=" + CodonConfiguration.workspace
						+ (isAccelerateSwiss() ? " --accelerate-swiss" : "")
						+ (fullBlast ? " --fullblast" : "")
						+ (optimizeIntegenicRegion? " --minimize-intergenic" : "")
						+ " -start=" + start
						+ " -end=" + end
						+ (!retreiveRNA ? " --ignore-rna" : "")
						+ (exportEmblFile != null ? " -export-emb=" + exportEmblFile.replaceAll (" ", "%20") : "")
						+ (exportCodonFile != null ? " -export-codon=" + exportCodonFile.replaceAll (" ", "%20") : "")
						+ " -overlaps-threehold=" + overlapThreeshold
						+ " -low-accuracy-threehold=" + accuracyFilter
						+ " -accuracy-tolerance=" + CodonConfiguration.tolerance
						+ " -accuracy-removing-overlap-during-blast=" + CodonConfiguration.accuracy_level_to_remove_overlap_during_blast
						+ " " + file);
			}
			else {
				String cmd =  "java -Xmx2548m "
						+ "-classpath lib/activation-1.1.jar:lib/animal-sniffer-annotations-1.14.jar:lib/antlr4-4.5.3.jar:lib/aopalliance-1.0.jar:lib/aopalliance-repackaged-2.4.0-b10.jar:lib/avro-1.7.6.jar:lib/avro-ipc-1.7.6.jar:lib/axis2-1.6.2.jar:lib/checker-compat-qual-2.0.0.jar:lib/classmate-0.8.0.jar:lib/commons-beanutils-1.9.2.jar:lib/commons-collections-3.2.1.jar:lib/commons-compress-1.12.jar:lib/commons-io-2.4.jar:lib/commons-lang-2.6.jar:lib/commons-lang3-3.6.jar:lib/commons-logging-1.2.jar:lib/commons-math3-3.6.1.jar:lib/commons-net-3.6.jar:lib/commons-text-1.8.jar:lib/error_prone_annotations-2.1.3.jar:lib/failsafe-1.0.3.jar:lib/FastInfoset-1.2.13.jar:lib/guava-24.1.1-jre.jar:lib/guice-4.0.jar:lib/guice-multibindings-4.0.jar:lib/hibernate-validator-5.0.1.Final.jar:lib/hk2-api-2.4.0-b10.jar:lib/hk2-locator-2.4.0-b10.jar:lib/hk2-utils-2.4.0-b10.jar:lib/httpclient-4.5.3.jar:lib/httpcore-4.4.6.jar:lib/httpmime-4.5.3.jar:lib/istack-commons-runtime-3.0.5.jar:lib/j2objc-annotations-1.1.jar:lib/jackson-annotations-2.9.5.jar:lib/jackson-core-2.9.5.jar:lib/jackson-core-asl-1.9.13.jar:lib/jackson-databind-2.9.4.jar:lib/jackson-mapper-asl-1.9.13.jar:lib/japi-1.0.30.jar:lib/javaparser-1.0.11.jar:lib/javassist-3.18.1-GA.jar:lib/javax.annotation-api-1.2.jar:lib/javax.inject-1.jar:lib/javax.inject-2.4.0-b10.jar:lib/javax.mail-1.6.1.jar:lib/javax.mail-api-1.6.1.jar:lib/javax.ws.rs-api-2.0.1.jar:lib/jaxb2-basics-1.11.1.jar:lib/jaxb2-basics-runtime-1.11.1.jar:lib/jaxb2-basics-tools-1.11.1.jar:lib/jaxb-api-2.3.0.jar:lib/jaxb-core-2.3.0.jar:lib/jaxb-impl-2.3.0.jar:lib/jaxb-runtime-2.3.0.jar:lib/jboss-logging-3.1.1.GA.jar:lib/jcl-over-slf4j-1.7.7.jar:lib/jcommander-1.32.jar:lib/jersey-client-2.17.jar:lib/jersey-common-2.17.jar:lib/jersey-guava-2.17.jar:lib/jersey-media-jaxb-2.17.jar:lib/jetty-6.1.26.jar:lib/jetty-util-6.1.26.jar:lib/jsr305-1.3.9.jar:lib/log4j-api-2.12.1.jar:lib/log4j-core-2.12.1.jar:lib/logback-classic-1.2.3.jar:lib/logback-core-1.2.3.jar:lib/metrics-core-3.2.2.jar:lib/netty-3.4.0.Final.jar:lib/noggit-0.8.jar:lib/osgi-resource-locator-1.0.1.jar:lib/paranamer-2.3.jar:lib/servlet-api-2.5-20081211.jar:lib/slf4j-api-1.7.25.jar:lib/snappy-java-1.0.5.jar:lib/solr-solrj-7.2.0.jar:lib/stax2-api-3.1.4.jar:lib/stax-ex-1.7.8.jar:lib/txw2-2.3.0.jar:lib/validation-api-1.1.0.Final.jar:lib/velocity-1.7.jar:lib/woodstox-core-asl-4.4.1.jar:lib/zookeeper-3.4.10.jar:./bin -Dswing.defaultlaf=javax.swing.plaf.nimbus.NimbusLookAndFeel"
						+ " codon.CommandLineBlasting"
						+ " -simultanous-blast=" + simultaneous
						+ " -workspace=" + CodonConfiguration.workspace
						+ (isAccelerateSwiss() ? " --accelerate-swiss" : "")
						+ (fullBlast ? " --fullblast" : "")
						+ (optimizeIntegenicRegion? " --minimize-intergenic" : "")
						+ " -start=" + start
						+ " -end=" + end
						+ (!retreiveRNA ? " --ignore-rna" : "")
						+ (exportEmblFile != null ? " -export-emb=" + exportEmblFile.replaceAll (" ", "%20") : "")
						+ (exportCodonFile != null ? " -export-codon=" + exportCodonFile.replaceAll (" ", "%20") : "")
						+ " -overlaps-threehold=" + overlapThreeshold
						+ " -low-accuracy-threehold=" + accuracyFilter
						+ " -accuracy-tolerance=" + CodonConfiguration.tolerance
						+ " -accuracy-removing-overlap-during-blast=" + CodonConfiguration.accuracy_level_to_remove_overlap_during_blast
						+ " " + file;
				proc = Runtime.getRuntime ().exec (cmd);
			}
		} catch (IOException e) {
			e.printStackTrace ();
		}
		ouP = new PrintStream (proc.getOutputStream ());
		new OutputControlerReader ().start ();
		new ErrorControlerReader ().start ();
	}
	
	Timer timer = new Timer (true);
	public BlastingMonitor (String file, int toProcess, int simultaneous, boolean accelerateSwiss, boolean fullBlast, boolean retreiveRNA, int start, int end) {
		if (!new File (file).exists ()) {
			System.err.println (file + " does not exist");
			return;
		}
		this.file = file;
		this.toProcess = this.toProcessAtStart = toProcess;
		this.simultaneous = simultaneous;
		this.setAccelerateSwiss (accelerateSwiss);
		this.fullBlast = fullBlast;
		this.start = start;
		this.end = end;
		this.retreiveRNA = retreiveRNA;
	}
	
	public void run () {
		lastHit = System.currentTimeMillis ();
		startProcess ();
		keyboardReader.start ();
		timer.schedule (new TimerTask () {
			public void run () {
				if (System.currentTimeMillis () - lastHit > 2500000) {
					System.err.println ("Process timout: restarting");
					lastHit = System.currentTimeMillis (); 
					timeoutRestart ++;
					kill ();
				}
			}
		}, 60000, 60000);
		try {
			proc.waitFor ();
			Thread.sleep (1000);
		} catch (InterruptedException e) {
			e.printStackTrace ();
		}
		while (!finalized) {
			System.err.println ("Process finalized before finalizing: restarting");
			startProcess ();
			deadRestart ++;
			fireRestarted ();
			try {
				Thread.sleep (1000);
				proc.waitFor ();
			} catch (InterruptedException e) {
				e.printStackTrace ();
			}
		}
	}
	
	public void setSimultaneousTries (int simultaneous) {
		this.simultaneous = simultaneous;
		ouP.println ("TRIES=" + simultaneous);
		ouP.flush ();
	}
	
	public void cancel () {
		kill ();
	}
	
	Thread keyboardReader = new Thread () {
		public void run () {
			BufferedReader reader = new BufferedReader (new InputStreamReader (System.in));
			String line = null;
			try {
				while ((line = reader.readLine ()) != null) {
					if (line.equals ("DIE")) {
						kill ();
						System.exit (0);
					}
					else if (line.equals ("RESTART")) {
						restartAll ();
					}
					else if (line.indexOf ("TRIES=") == 0) {
						simultaneous = Integer.parseInt (line.substring (6));
						setSimultaneousTries (simultaneous);
					}
					else {
						System.out.println ("Blasting: " + blasting + " / " + simultaneous);
						System.out.println ("Ended: " + ended);
						System.out.println ("Remaining: " + toProcess);
						System.out.println ("Previously done: " + previouslyDone);
						System.out.println ("Ignored: " + ignored);
						System.out.println ("Max memory load: " + maxLoad);
						System.out.println ("Dead restarts: " + deadRestart);
						System.out.println ("Overload memory restarts: " + memoryRestart);
						System.out.println ("Timeout restarts: " + timeoutRestart);
						System.out.println (status);
					}
				}
			} catch (IOException e) {
				e.printStackTrace ();
			}
		};
	};
	
	static String exportEmblFile = null;
	static String exportCodonFile = null;
	static int accuracyFilter = -1;
	static int overlapThreeshold = -1;
	static boolean optimizeIntegenicRegion = false;
	static boolean verbose = false;
	public static boolean blastErrors = false;
	
	public static void main (String [] args) {
		/*if (args.length == 0) {
			new BlastingMonitor ("./workspace/M26365.gb", -1, 10, false, false, true, -1, -1).start ();
			return;
		}*/
		standaloneApplication = true;
		CodonConfiguration.ui = false;
		CodonConfiguration.workspace = "./workspace";
		int threadNumber = 1;
		boolean fullBlast = false;
		boolean accelerateWithSwiss = false;
		boolean retreiveRNA = true;
		int start = -1;
		int end = -1;
		
		for (int i = 0; i < args.length; i++) {
			if (args [i].equalsIgnoreCase ("--verbose")) verbose = true;
			if (args [i].equalsIgnoreCase ("--see-blast-error")) blastErrors = true;
			if (args [i].equalsIgnoreCase ("--fullblast")) fullBlast = true;
			if (args [i].equalsIgnoreCase ("--ignore-rna")) retreiveRNA = false;
			if (args [i].equalsIgnoreCase ("--accelerate-swiss")) accelerateWithSwiss = true;
			if (args [i].indexOf ("-simultanous-blast=") == 0) threadNumber = CommandLineBlasting.safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-start=") == 0) start = CommandLineBlasting.safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-end=") == 0) end = CommandLineBlasting.safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-workspace=") == 0) CodonConfiguration.workspace = args [i].substring (args [i].indexOf ("=") + 1);
			if (args [i].indexOf ("-export-emb=") == 0) exportEmblFile = args [i].substring (args [i].indexOf ("=") + 1);
			if (args [i].indexOf ("-export-codon=") == 0) exportCodonFile = args [i].substring (args [i].indexOf ("=") + 1);
			if (args [i].indexOf ("-overlapss-threehold") == 0) overlapThreeshold = CommandLineBlasting.safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-low-accuracy-threehold=") == 0) accuracyFilter = CommandLineBlasting.safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-accuracy-tolerance=") == 0) CodonConfiguration.tolerance = CommandLineBlasting.safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("-accuracy-removing-overlap-during-blast=") == 0) CodonConfiguration.accuracy_level_to_remove_overlap_during_blast = CommandLineBlasting.safeParamToIntConversion (args [i]);
			if (args [i].indexOf ("--minimize-intergenic") == 0) optimizeIntegenicRegion = true;
			if (args [i].indexOf ("-status-file=") == 0) statusFile = args [i].substring (args [i].indexOf ("=") + 1);
		}
		
		if (args.length == 0) {
			System.err.println ("The command line must end by an existing fasta file or embl file");
			System.exit (0);
		}
		
		String file = args [args.length - 1];
		if (!new File (file).exists ()) {
			System.err.println ("The file " + file + " does not exist");
			System.exit (0);
		}
		new BlastingMonitor (file, -1, threadNumber, accelerateWithSwiss, fullBlast, retreiveRNA, start, end).start ();
	}

	public boolean isAccelerateSwiss () {
		return accelerateSwiss;
	}

	public void setAccelerateSwiss (boolean accelerateSwiss) {
		this.accelerateSwiss = accelerateSwiss;
		if (ouP == null) return;
		ouP.println ("SWISS=" + accelerateSwiss);
		ouP.flush ();
	}
}
