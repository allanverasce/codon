package codon.recognizer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import codon.data.Codon;
import codon.data.CodonConfiguration;
import codon.data.OrganismFamilly;
import codon.data.GeneFeature;
import codon.data.LexicographicTree;
import codon.data.ORF;
import codon.data.Organism;
import codon.data.Overlap;
import codon.data.ProductProtein;
import codon.fsm.FSMOrfRetriever;
import splash.WaitingSplash;
import splash.WaitingSplashThread;

public class OrfRecognizer {
	
	public final static int V = 0;
	public final static int A = 1;
	public final static int T = 2;
	public final static int C = 3;
	public final static int AT = 4;
	public final static int TA = 5;
	public final static int TG = 6;
	public final static int TT = 7;
	public final static int TC = 8;
	public final static int CT = 9;
	public final static int CA = 10;
	public final static int G = 11;
	public final static int AA = 12;
	public final static int GT = 13;
	public final static int GA = 14;
	
	static Hashtable <String, String> acidoBaseMap = new Hashtable <String, String> ();
	public static Hashtable <String, String> baseAcidoMap = new Hashtable <String, String> ();
	
	int matrix [][] = {
			{V, 'A', A, NONE},
			{V, 'G', G, NONE},
			{V, 'T', T, NONE},
			{V, 'C', C, NONE},
			
			{A, 'A', AA, NONE},
			{A, 'G', G, NONE},
			{A, 'T', AT, NONE},
			{A, 'C', C, NONE},
			
			{T, 'A', TA, NONE},
			{T, 'G', TG, NONE},
			{T, 'T', TT, NONE},
			{T, 'C', TC, NONE},
			
			{C, 'A', CA, NONE},
			{C, 'G', G, NONE},
			{C, 'T', CT, NONE},
			{C, 'C', C, NONE},
			
			{AT, 'A', TA, NONE},
			{AT, 'G', TG, START},
			{AT, 'T', TT, NONE},
			{AT, 'C', TC, NONE},
			
			{TA, 'A', AA, STOP},
			{TA, 'G', G, STOP},
			{TA, 'T', AT, NONE},
			{TA, 'C', C, VALINA_R},
			
			{TG, 'A', GA, STOP},
			{TG, 'G', G, NONE},
			{TG, 'T', GT, NONE},
			{TG, 'C', C, NONE},
			
			{TT, 'A', TA, STOP_R},
			{TT, 'G', TG, LEUCINA},
			{TT, 'T', TT, NONE},
			{TT, 'C', TC, NONE},
			
			{TC, 'A', CA, STOP_R},
			{TC, 'G', G, NONE},
			{TC, 'T', CT, NONE},
			{TC, 'C', C, NONE},
			
			{CT, 'A', TA, STOP_R},
			{CT, 'G', TG, LEUCINA},
			{CT, 'T', TT, LEUCINA},
			{CT, 'C', TC, LEUCINA},
			
			{CA, 'A', AA, LEUCINA_R},
			{CA, 'G', G, LEUCINA_R},
			{CA, 'T', AT, START_R},
			{CA, 'C', C, VALINA_R},
			
			{G, 'A', GA, NONE},
			{G, 'G', G, NONE},
			{G, 'T', GT, NONE},
			{G, 'C', C, NONE},
			
			{AA, 'A', AA, NONE},
			{AA, 'G', G, LEUCINA_R},
			{AA, 'T', AT, NONE},
			{AA, 'C', C, VALINA_R},
			
			{GT, 'A', TA, VALINA},
			{GT, 'G', TG, VALINA},
			{GT, 'T', TT, VALINA},
			{GT, 'C', TC, VALINA},
			
			{GA, 'A', AA, NONE},
			{GA, 'G', G, LEUCINA_R},
			{GA, 'T', AT, NONE},
			{GA, 'C', C, VALINA_R},
	};
	
	static {
		baseAcidoMap.put ("TTT", "F");
		baseAcidoMap.put ("TTC", "F");
		//acidoBaseMap.put ("F", "TT(T|C)");
		acidoBaseMap.put ("F", "TT.");
		baseAcidoMap.put ("TTA", "L");
		baseAcidoMap.put ("TTG", "L");		
				
		baseAcidoMap.put ("CTT", "L");
		baseAcidoMap.put ("CTC", "L");
		baseAcidoMap.put ("CTA", "L");
		baseAcidoMap.put ("CTG", "L");
		//acidoBaseMap.put ("L", "(CT.|TT(A|G))");
		acidoBaseMap.put ("L", ".T.");
		
		baseAcidoMap.put ("ATT", "I");
		baseAcidoMap.put ("ATC", "I");
		baseAcidoMap.put ("ATA", "I");
		//acidoBaseMap.put ("I", "AT(T|C|A)");
		acidoBaseMap.put ("I", "AT.");
		baseAcidoMap.put ("ATG", "M");
		acidoBaseMap.put ("M", "ATG");
		
		baseAcidoMap.put ("GTT", "V");
		baseAcidoMap.put ("GTC", "V");
		baseAcidoMap.put ("GTA", "V");
		baseAcidoMap.put ("GTG", "V");
		acidoBaseMap.put ("V", "GT.");
		
		///
		
		baseAcidoMap.put ("TCT", "S");
		baseAcidoMap.put ("TCC", "S");
		baseAcidoMap.put ("TCA", "S");
		baseAcidoMap.put ("TCG", "S");
//		acidoBaseMap.put ("S", "(AG(T|C)|TC.)");
		acidoBaseMap.put ("S", "...");
		
		baseAcidoMap.put ("CCT", "P");
		baseAcidoMap.put ("CCC", "P");
		baseAcidoMap.put ("CCA", "P");
		baseAcidoMap.put ("CCG", "P");
		acidoBaseMap.put ("P", "CC.");
		
		baseAcidoMap.put ("ACT", "T");
		baseAcidoMap.put ("ACC", "T");
		baseAcidoMap.put ("ACA", "T");
		baseAcidoMap.put ("ACG", "T");
		acidoBaseMap.put ("T", "AC.");
		
		baseAcidoMap.put ("GCT", "A");
		baseAcidoMap.put ("GCC", "A");
		baseAcidoMap.put ("GCA", "A");
		baseAcidoMap.put ("GCG", "A");
		acidoBaseMap.put ("A", "GC.");
		
		///
		
		baseAcidoMap.put ("TAT", "Y");
		baseAcidoMap.put ("TAC", "Y");
//		acidoBaseMap.put ("Y", "TA(T|C)");
		acidoBaseMap.put ("Y", "TA.");
		baseAcidoMap.put ("TAA", "#");
		baseAcidoMap.put ("TAG", "+");
		
		baseAcidoMap.put ("CAT", "H");
		baseAcidoMap.put ("CAC", "H");
//		acidoBaseMap.put ("H", "CA(T|C)");
		acidoBaseMap.put ("H", "CA.");
		baseAcidoMap.put ("CAA", "Q");
		baseAcidoMap.put ("CAG", "Q");
//		acidoBaseMap.put ("Q", "CA(A|G)");
		acidoBaseMap.put ("Q", "CA.");
		
		
		baseAcidoMap.put ("AAT", "N");
		baseAcidoMap.put ("AAC", "N");
//		acidoBaseMap.put ("N", "AA(T|C)");
		acidoBaseMap.put ("N", "AA.");
		baseAcidoMap.put ("AAA", "K");
		baseAcidoMap.put ("AAG", "K");
//		acidoBaseMap.put ("K", "AA(A|G)");
		acidoBaseMap.put ("K", "AA.");
		
		baseAcidoMap.put ("GAT", "D");
		baseAcidoMap.put ("GAC", "D");
//		acidoBaseMap.put ("D", "GA(T|C)");
		acidoBaseMap.put ("D", "GA.");
		baseAcidoMap.put ("GAA", "E");
		baseAcidoMap.put ("GAG", "E");
//		acidoBaseMap.put ("E", "AA(A|G)");
		acidoBaseMap.put ("E", "AA.");
		
		///
		
		baseAcidoMap.put ("TGT", "C");
		baseAcidoMap.put ("TGC", "C");
//		acidoBaseMap.put ("Y", "TG(T|C)");
		acidoBaseMap.put ("C", "TG.");
		baseAcidoMap.put ("TGA", "*");
		baseAcidoMap.put ("TGG", "W");
		acidoBaseMap.put ("W", "TGG");
		
		baseAcidoMap.put ("CGT", "R");
		baseAcidoMap.put ("CGC", "R");
		baseAcidoMap.put ("CGA", "R");
		baseAcidoMap.put ("CGG", "R");
//		acidoBaseMap.put ("R", "(CG.|AG(A|G))");
		acidoBaseMap.put ("R", ".G.");
		
		baseAcidoMap.put ("AGT", "S");
		baseAcidoMap.put ("AGC", "S");
		baseAcidoMap.put ("AGA", "R");
		baseAcidoMap.put ("AGG", "R");
		
		baseAcidoMap.put ("GGT", "G");
		baseAcidoMap.put ("GGC", "G");
		baseAcidoMap.put ("GGA", "G");
		baseAcidoMap.put ("GGG", "G");
		acidoBaseMap.put ("G", "GG.");
		
		acidoBaseMap.put ("X", "...");
		acidoBaseMap.put ("U", "...");
	}
	
	
	public int numberA = 0;
	public int numberT = 0;
	public int numberC = 0;
	public int numberG = 0;
	
	public final static int GC_WINDOW_SIZE = 192;
	
	public final static int NONE = -1;
	public final static int START = 0;
	public final static int STOP = 1;
	public final static int VALINA = 2;
	public final static int LEUCINA = 2;
        
	public final static int START_R = 3;
	public final static int STOP_R = 4;
    public final static int VALINA_R = 5;
    public final static int LEUCINA_R = 5;
    
    public LexicographicTree duplicatedEntries = null;
	public LexicographicTree duplicatedGenes = null;
        
	public String dataDestDirS = null;
	public String dataSrcFileS = null;
	String codonsType [] = {"START", "STOP", "START_R", "STOP_R"};
	
	@SuppressWarnings ("unchecked")
	public
	ArrayList <Codon> [] startStopPorFrame =  new ArrayList [6] ;
	@SuppressWarnings ("unchecked")
	public ArrayList <ORF> [] orfsPorFrame = new ArrayList [6];
	public List <ORF> rnaOrfs = new ArrayList <> ();
//	public List <ORF> tRNAOrfs = new ArrayList <> ();
	
	ArrayList <OrfRecognizerListener> listeners = new ArrayList <> ();
	public void addOrfRecognizerListener (OrfRecognizerListener orfRecognizerListener) {
		listeners.add (orfRecognizerListener);
	}
	public void removeOrfRecognizerListener (OrfRecognizerListener orfRecognizerListener) {
		listeners.remove (orfRecognizerListener);
	}
	public void fireSuperpositionsAndFreeSpansUpdate () {
		for (OrfRecognizerListener l: listeners) l.superpositionsAndFreeSpansUpdate ();
	}
	
	public OrfRecognizer () {
		clear ();
	}
	
	public void addChar (char c) {
		switch (c) {
			case 'a': c = 'A'; break;
			case 't': c = 'T'; break;
			case 'g': c = 'G'; break;
			case 'c': c = 'C'; break;
			default: break;
		}
		int trans = currentState * 4;
		switch (c) {
		case 'A':
			numberA ++;
			break;
		case 'G':
			numberG ++;
			trans += 1;
			break;
		case 'T':
			numberT ++;
			trans += 2;
			break;
		case 'C':
			numberC ++;
			trans += 3;
			break;
		default:
			break;
		}
		
		int [] transData = matrix [trans];
		currentState = transData [2];
		if (transData [3] >= 0) addIntoFrame (index - 2, transData [3]);
		index ++;
	}
	
	private int currentState = 0;
	private int index = 0;
	
	public void loadFastaString (String destRep, String bases) {
		File f = new File (destRep);
		currentState = 0;
		index = 0;
		this.dataSrcFileS = destRep;
		this.dataDestDirS = f.getName ();
		File dir = new File (CodonConfiguration.workspace + dataDestDirS + "/orfs");
		if (!dir.exists ()) dir.mkdirs ();
		for (int i = 0; i < bases.length (); i++) {
			addChar (bases.charAt (i));
		}
		wholeSequence = new StringBuffer (bases);
		prepare ();
	}
	
	public void load (String file) {
		File f = new File (file);
		currentState = 0;
		index = 0;
		this.dataSrcFileS = file;
		this.dataDestDirS = f.getName ();
		this.dataDestDirS = dataDestDirS.substring (0, dataDestDirS.lastIndexOf ('.')); 
		File dir = new File (CodonConfiguration.workspace + dataDestDirS + "/orfs");
		if (!dir.exists ()) dir.mkdirs ();
		try {
			BufferedReader reader = new BufferedReader (new FileReader (dataSrcFileS));
			String l = reader.readLine ();
			while ((l = reader.readLine ()) != null) {
				for (int i = 0; i < l.length (); i++) {
					addChar (l.charAt (i));
				}
			}
			reader.close ();
		} catch (FileNotFoundException e) {
			e.printStackTrace ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
		
		prepare ();
		loadBaseSequence ();
	}
	
	public void clear () {
		numberA = 0;
		numberT = 0;
		numberC = 0;
		numberG = 0;
		for (int i = 0; i < metiolinaStart.length; i++) metiolinaStart [i] = -1;
		duplicatedEntries = new LexicographicTree ();
		duplicatedGenes = new LexicographicTree ();
		for (int frame = 0; frame < 6; frame++) {
			startStopPorFrame [frame] = new ArrayList <Codon> ();
			orfsPorFrame [frame] = new ArrayList <ORF> ();
		}
		rnaOrfs.clear ();
//		tRNAOrfs.clear ();
	}
	
	private Codon lastSSCodon [] = new Codon [6];
	

	int metiolinaStart [] = new int [6];
	public void addIntoFrame (int index, int type) {
		int frameIndex = (index  + 1) % 3 + 3 * (type / 3);
		ArrayList <Codon> startStops = startStopPorFrame [frameIndex];
		List <ORF> orfs = orfsPorFrame [frameIndex];
		
		if (type < 3) {
			if ((type == START || type == VALINA) && (orfs.size () == 0 || lastSSCodon [frameIndex].type == STOP)) {
				if (orfs.size () != 0 && orfs.get (orfs.size () - 1).length < CodonConfiguration.min_orf_length) orfs.remove (orfs.size () - 1);
				orfs.add (new ORF (this, index, -1, frameIndex));
				if (type == START) metiolinaStart [frameIndex] = index;
			}
			if (metiolinaStart [frameIndex] == -1 && type == START) {
				metiolinaStart [frameIndex] = index;
			}
			if (type == STOP && orfs.size () > 0 && (lastSSCodon [frameIndex].type == START || lastSSCodon [frameIndex].type == VALINA)) {
				if (metiolinaStart [frameIndex] != - 1 && index - metiolinaStart [frameIndex] > CodonConfiguration.min_orf_length) orfs.get (orfs.size () - 1).setOriginalStart (metiolinaStart [frameIndex]);
				orfs.get (orfs.size () - 1).setOriginalStop (index);
				metiolinaStart [frameIndex] = -1;
			}
		}
		else {
			if (lastSSCodon [frameIndex] != null && lastSSCodon [frameIndex].type == STOP_R && (type == START_R || type == VALINA_R)) {
				if (orfs.size () != 0 && orfs.get (orfs.size () - 1).length < CodonConfiguration.min_orf_length) orfs.remove (orfs.size () - 1);
				orfs.add (new ORF (this, lastSSCodon [frameIndex].index, index, frameIndex));
				metiolinaStart [frameIndex] = -1;
			}
			else {
				if (orfs.size () > 0 && type == START_R) {
					orfs.get (orfs.size () - 1).setOriginalStop (index);
					metiolinaStart [frameIndex] = index;
				}
				if (orfs.size () > 0 && type == VALINA_R && (orfs.get (orfs.size () - 1).length < CodonConfiguration.min_orf_length || metiolinaStart [frameIndex] == -1)) orfs.get (orfs.size () - 1).setOriginalStop (index);
			}
		}
    	if (frameIndex < 3 && (type == START || type == STOP || type == VALINA) || frameIndex > 2 && (type == START_R || type == STOP_R || type == VALINA_R)) {
    		lastSSCodon [frameIndex] = new Codon (type, index);
    		startStops.add (lastSSCodon [frameIndex]);
    	}
	}
	
	private void prepare () {
		for (int i = 0; i < 3; i++) {
			if (orfsPorFrame [i].get (orfsPorFrame [i].size () - 1).stop == -1) orfsPorFrame [i].remove (orfsPorFrame [i].size () - 1);
		}
		for (int i = 0; i < orfsPorFrame.length; i++) {
			List <ORF> orfs = orfsPorFrame [i];
			for (ORF orf: orfs) {
				orf.setRemoved (orf.length < CodonConfiguration.min_orf_length);
				orf.rename ();
			}
		}
	}
	
	public void reboundStarts () {
		for (List <ORF> orfs: orfsPorFrame) {
			for (ORF orf: orfs) {
				if (orf.isRemoved ()) continue;
				orf.autoRebound ();
			}
		}
	}
	
	public void selectBetterEntriesThatFit () {
		total = 0;
		total = 0;
		for (List <ORF> orfs: orfsPorFrame) {
			for (ORF orf: orfs) {
				if (!orf.isRemoved ()) total ++;				
			}
		}
		WaitingSplashThread th = new WaitingSplashThread () {
			public void run () {
				int count = 0;
				for (List <ORF> orfs: orfsPorFrame) {
					for (ORF orf: orfs) {
						if (splash != null) splash.changeStatus (count + " / " + total, count);
						if (mustStop) return;
						if (orf.isRemoved ()) continue;
						if (orf.isGene () && orf.getAccuracy () == 100) continue;
						orf.selectBetterEntryThatFit ();
					}
				}
				if (splash != null) splash.setVisible (false);
			}
		};
		if (CodonConfiguration.ui) {
			splash = new WaitingSplash (null, "Selecting the entries that better fit", "0 / " + total, total, th, true);
			th.start ();
			splash.setVisible (true);
		}
		else {
			th.run ();
		}
	}
	
	public void selectBetterEntriesMaximizingIntergenic () {
		total = 0;
		for (List <ORF> orfs: orfsPorFrame) {
			for (ORF orf: orfs) {
				if (!orf.isRemoved ()) total ++;				
			}
		}
		WaitingSplashThread th = new WaitingSplashThread () {
			public void run () {
				int count = 0;
				for (List <ORF> orfs: orfsPorFrame) {
					for (ORF orf: orfs) {
						if (splash != null) splash.changeStatus (count + " / " + total, count);
						if (mustStop) return;
						if (orf.isRemoved ()) continue;
						orf.selectBetterEntryMaximizingIntergenicRegion ();
					}
				}
				if (splash != null) splash.setVisible (false);
			}
		};
		if (CodonConfiguration.ui) {
			splash = new WaitingSplash (null, "Selecting the entries that better fit", "0 / " + total, total, th, true);
			th.start ();
			splash.setVisible (true);
		}
		else {
			th.run ();
		}
	}
	
	private int total = 0;
	public void selectBetterEntriesFulfillIntergenic (int tolerance) {
		total = 0;
		total = 0;
		for (List <ORF> orfs: orfsPorFrame) {
			for (ORF orf: orfs) {
				if (!orf.isRemoved ()) total ++;				
			}
		}
		WaitingSplashThread th = new WaitingSplashThread () {
			public void run () {
				int count = 0;
				for (List <ORF> orfs: orfsPorFrame) {
					for (ORF orf: orfs) {
						if (splash != null) splash.changeStatus (count + " / " + total + " (" + orf.id + ")", count);
						if (mustStop) return;
						if (orf.isRemoved ()) continue;
						count ++;
						if (orf.isForward ()) {
							int ind = getCodonIndex (intergenicRegions, orf.start - 1);
							int startSpan = intergenicRegions.get (ind).index;
							int stopSpan = intergenicRegions.get (ind + 1).index;
							if (stopSpan == orf.start) {
								orf.selectBetterEntryFulfill (orf.length + stopSpan - startSpan + tolerance);
							}
						}
						else {
							int ind = getCodonIndex (intergenicRegions, orf.stop + 1);
							int startSpan = intergenicRegions.get (ind).index;
							int stopSpan = intergenicRegions.get (ind + 1).index;
							if (startSpan == orf.stop) {
								orf.selectBetterEntryFulfill (orf.length + stopSpan - startSpan + tolerance);
							}
		//					int ind = getCodonIndex (freeSpans, orf.stop + 1);
		//					orf.selectBetterEntryFulfill (orf.lenght + freeSpans.get (ind).index - freeSpans.get (ind + 1).index + tolerance);
						}
					}
				}
				if (splash != null) splash.setVisible (false);
			}
		};
		if (CodonConfiguration.ui) {
			splash = new WaitingSplash (null, "Selecting the entries that better fulfill the integenic regions", "0 / " + total, total, th, true);
			th.start ();
			splash.setVisible (true);
		}
		else {
			th.run ();
		}
	}
	
//	private void removeIntersectionInOrf0 (ORF orf0, ORF orf1, int intersection) {
//		if (intersection > orf0.lenght / 2) orf0.setRemoved (true);
//		else if (!orf0.isReverse () && !orf1.isReverse ()) {
//			if (orf0.start > orf1.start) orf0.castrating (intersection);
//			//else orf0.beheading (intersection);
//		}
//		else if (orf0.isReverse () && orf1.isReverse ()) {
//			if (orf0.start < orf1.start) orf0.castrating (intersection);
//			//else orf0.beheading (intersection);
//		}
//		else {
//			//if (!orf0.reverse && orf0.start < orf1.start) orf0.beheading (intersection);
//			if (!orf0.isReverse () && orf0.start > orf1.start) orf0.castrating (intersection);
//			if (orf0.isReverse () && orf0.start < orf1.start) orf0.castrating (intersection);
//			//if (orf0.reverse && orf0.start > orf1.start) orf0.beheading (intersection);
//		}
//	}
	
	public ArrayList <Codon> intergenicRegions = new ArrayList <> ();
	public ArrayList <Codon> overlapSpans = new ArrayList <> ();
	ArrayList <Overlap> overlaps = new ArrayList <> ();
	
	@SuppressWarnings ("unchecked")
	public ArrayList <Codon> [] computeOrfFreeAndOrfSuperpositionSpans (int overlapTolerance) {
		int currentIndex [] = new int [] {0, 0, 0, 0, 0, 0};
		ORF currentOrfs [] = new ORF [] {
			orfsPorFrame [0].size () == 0 ? null : orfsPorFrame [0].get (0), 
			orfsPorFrame [1].size () == 0 ? null : orfsPorFrame [1].get (0), 
			orfsPorFrame [2].size () == 0 ? null : orfsPorFrame [2].get (0),
			orfsPorFrame [3].size () == 0 ? null : orfsPorFrame [3].get (0),
			orfsPorFrame [4].size () == 0 ? null : orfsPorFrame [4].get (0),
			orfsPorFrame [5].size () == 0 ? null : orfsPorFrame [5].get (0),
		};
		
		intergenicRegions.clear ();
		overlapSpans.clear ();
		overlaps.clear ();
		intergenicRegions.add (new Codon (START, 0));
		
		for (int i = 0; i < orfsPorFrame.length; i++) {
			List <ORF> orfs = orfsPorFrame [i];
			for (int j = 0; j < orfs.size (); j++) {
				orfs.get (j).overlaps = 0;
			}
		}
		
		int lastStart = Integer.MAX_VALUE;
		int lastStop = 0;
		int lastSobrpositionStop = -1;
		ORF lastOrf = null;
		
		int total = 0;
		for (int i = 0; i < 6; i++) total += orfsPorFrame [i].size ();
		
		while (total > 0) {
			int min = -1;
			int stMin = Integer.MAX_VALUE;
			for (int i = 0; i < currentOrfs.length; i++) {
				if (currentIndex [i] < orfsPorFrame [i].size () && currentOrfs [i].start < stMin) {
					stMin = currentOrfs [i].start;
					min = i;
				}
			}
			ORF orf = currentOrfs [min];
			
			if (!orf.isRemoved ()) {
				if (orf.start > lastStart && orf.start < lastStop) {
					if (overlapSpans.size () == 0 || orf.start > lastSobrpositionStop) {
						if (overlapSpans.size () != 0) {
							overlapSpans.add (new Codon (STOP, lastSobrpositionStop));
							overlaps.get (overlaps.size () - 1).end = lastSobrpositionStop;
						}
						overlapSpans.add (new Codon (START, orf.start));
						lastSobrpositionStop = Math.min (lastStop, orf.stop);
						Overlap ov = new Overlap (orf.start, lastSobrpositionStop);
						ov.addOrf (lastOrf);
						lastOrf.overlaps ++;
						ov.addOrf (orf);
						orf.overlaps ++;
						overlaps.add (ov);
						lastStop = Math.max (lastStop, orf.stop);
					}
					else { 
						if (orf.stop < lastStop && orf.stop > lastSobrpositionStop) {
							lastSobrpositionStop = orf.stop;
							overlaps.get (overlaps.size () - 1).addOrf (orf);
						}
						else if (orf.stop > lastStop) {
							lastSobrpositionStop = lastStop;
							lastStop = orf.stop;
							overlaps.get (overlaps.size () - 1).addOrf (orf);
						}
					}
				}
				
				if (intergenicRegions.size () == 0 || orf.start > intergenicRegions.get (intergenicRegions.size () - 1).index) {
//					if (orf.start - freeSpans.get (freeSpans.size () - 1).index > ORF.MIN_LENGHT - 20) {
						intergenicRegions.add (new Codon (STOP, orf.start));
						lastStart = orf.start;
						intergenicRegions.add (new Codon (START, orf.stop));
						lastStop = orf.stop;
//					}
//					else {
//						freeSpans.get (freeSpans.size () - 1).index = orf.stop;
//						lastStop = orf.stop;
//					}
				}
				else if (orf.stop > intergenicRegions.get (intergenicRegions.size () - 1).index) {
					intergenicRegions.get (intergenicRegions.size () - 1).index = orf.stop;
					lastStop = orf.stop;
				}
				
				lastOrf = orf;
			}
			
			currentIndex [min] ++;
			if (currentIndex [min] < orfsPorFrame [min].size ()) currentOrfs [min] =  orfsPorFrame [min].get (currentIndex [min]);
			total --;
		}
		intergenicRegions.add (new Codon (STOP, wholeSequence.length ()));
		if (overlapSpans.size () != 0) {
			overlapSpans.add (new Codon (STOP, lastSobrpositionStop));
			overlaps.get (overlaps.size () - 1).end = lastSobrpositionStop;
		}
//		for (Codon cod: intergenicRegions) {
//			if (cod.type == START) cod.index -= 1;
//			if (cod.type == STOP) cod.index += 1;
//		}
		ArrayList <Codon> spans = new ArrayList <Codon> ();
		for (int i = 0; i < overlapSpans.size (); i += 2) {
			Codon start = overlapSpans.get (i);
			Codon stop = overlapSpans.get (i + 1);
			if (stop.index - start.index > overlapTolerance) {
				spans.add (start);
				spans.add (stop);
			}
		}
		overlapSpans = spans;
		fireSuperpositionsAndFreeSpansUpdate ();
		return new ArrayList [] {intergenicRegions, overlapSpans};
	}
	
	public ArrayList <ORF> getOverlapedOrfs () {
		ArrayList<ORF> overlapedOrfs = new ArrayList <> ();
		int currentIndex [] = new int [] {0, 0, 0, 0, 0, 0};
		ORF currentOrfs [] = new ORF [] {orfsPorFrame [0].get (0), orfsPorFrame [1].get (0), orfsPorFrame [2].get (0), 
										 orfsPorFrame [3].get (0), orfsPorFrame [4].get (0), orfsPorFrame [5].get (0)};
		
		int total = 0;
		for (int i = 0; i < 6; i++) total += orfsPorFrame [i].size ();
		ORF lastOrf = null;
		boolean added = false;
		
		while (total > 0) {
			int min = -1;
			int stMin = Integer.MAX_VALUE;
			for (int i = 0; i < currentOrfs.length; i++) {
				if (currentIndex [i] < orfsPorFrame [i].size () && currentOrfs [i].start < stMin) {
					stMin = currentOrfs [i].start;
					min = i;
				}
			}
			ORF orf = currentOrfs [min];
			if (!orf.isRemoved ()) {
				if (lastOrf != null) {
					if (orf.start > lastOrf.start && orf.start < lastOrf.stop) {
						if (!added) overlapedOrfs.add (lastOrf);
						overlapedOrfs.add (orf);
						added = true;
					}
					else {
						added = false;
					}
				}
				lastOrf = orf;
			}
			
			currentIndex [min] ++;
			if (currentIndex [min] < orfsPorFrame [min].size ()) currentOrfs [min] =  orfsPorFrame [min].get (currentIndex [min]);
			total --;
		}

		for (Codon cod: intergenicRegions) {
			if (cod.type == START) cod.index -= 1;
			if (cod.type == STOP) cod.index += 1;
		}
		return overlapedOrfs;
	}
	
	public int getCodonIndex (List <Codon> codons, int index) {
		return getCodonIndex (codons, 0, codons.size () - 1, index);
	}
	
	private int getCodonIndex (List <Codon> codons, int start, int stop, int index) {
		if (stop - start <= 1) return start;
		int mid = (stop + start) / 2;
		Codon codon = codons.get (mid);
		if (codon.index == index) return mid;
		if (codon.index > index) return getCodonIndex (codons, start, mid, index);
		else return getCodonIndex (codons, mid, stop, index);
	}
	
	public void clearRemoved () {
		for (int i = 0; i < 6; i++) for (ORF orf: orfsPorFrame [i]) orf.setRemoved (orf.length < CodonConfiguration.min_orf_length);
	}
	
//	public void removeTinyUnspecifyingSequence () {
//		for (int i = 0; i < 6; i++) {
//			for (ORF orf: orfsPorFrame [i]) {
//				if (orf.lenght < 120 && (orf.getSpecificity () < 90 || orf.isUncharacterized () && !orf.isGene ())) orf.setRemoved (true);
//			}
//		}
//	}
	
	public int candidateOrfCount () {
		int total = 0;
		for (int i = 0; i < 6; i++) for (ORF orf: orfsPorFrame [i]) if (!orf.isRemoved ()) total ++;
		return total;
	}
	
	private void removeSuperpositions (List <ORF> orfs0, List <ORF> orfs1) {
		Iterator <ORF> it1 = orfs1.iterator ();
		if (!it1.hasNext ()) return;
		ORF orf1 = it1.next ();
		while (it1.hasNext () && (orf1.isRemoved () || !orf1.isBlasted ())) {
			orf1 = it1.next ();
		}
		for (int i = 0; i < orfs0.size (); i ++) {
			if (orf1.isRemoved ()) return;
			ORF orf0 = orfs0.get (i);
			if (orf0.isRemoved ()) continue;
			while (true) {
				if (orf1.start > orf0.stop) break;
				int intersection = orf0.intersection (orf1);
				if (intersection > orf0.length * 30 / 100) {
					if (!orf0.isGene () && orf1.isGene () && orf1.getAccuracy () >= orf0.getAccuracy () - CodonConfiguration.tolerance) orf0.setRemoved (true);
					else if (!orf1.isGene () && orf0.isGene () && orf0.getAccuracy () >= orf1.getAccuracy () - CodonConfiguration.tolerance) orf1.setRemoved (true);
					else if (orf0.isGene () && orf1.isGene () && !orf0.sizesFit () && orf1.sizesFit ())  orf0.setRemoved (true);
					else if (orf1.isGene () && orf0.isGene () && !orf1.sizesFit () && orf0.sizesFit ())  orf1.setRemoved (true);
					else if (orf0.isUncharacterized () && !orf1.isUncharacterized ()) orf0.setRemoved (true);
					else if (orf1.isUncharacterized () && !orf0.isUncharacterized ()) orf1.setRemoved (true);
					else if (orf0.isUncharacterized () && orf1.isUncharacterized ()) {
						if (!orf0.sizesFit () && orf1.sizesFit ()) orf0.setRemoved (true);
						else if (!orf1.sizesFit () && orf0.sizesFit ()) orf1.setRemoved (true);
						else if (orf0.inside (orf1)) orf0.setRemoved (true);
						else if (orf1.inside (orf0)) orf1.setRemoved (true);
						else if (orf0.overlaps > orf1.overlaps) orf0.setRemoved (true);
					}
					else if (!orf0.isUncharacterized () && !orf1.isUncharacterized ()) {
						if (!orf0.sizesFit () && orf1.sizesFit ()) orf0.setRemoved (true);
						else if (!orf1.sizesFit () && orf0.sizesFit ()) orf1.setRemoved (true);
						else if (orf0.inside (orf1)) orf0.setRemoved (true);
						else if (orf1.inside (orf0)) orf1.setRemoved (true);
						else if (orf0.overlaps > orf1.overlaps) orf0.setRemoved (true);
					}
				}
				if (orf0.isRemoved () || intersection > 0 && orf0.stop < orf1.stop) break;
				if (it1.hasNext ()) orf1 = it1.next ();
				while (it1.hasNext () && (orf1.isRemoved () || !orf1.isBlasted ())) {
					orf1 = it1.next ();
				}
				if (orf1.isRemoved () || !orf1.isBlasted () || !it1.hasNext ()) break;
			}
		}
	}
	
	WaitingSplash splash;
	public void removeOverlaps (int tolerance) {
		computeOrfFreeAndOrfSuperpositionSpans (tolerance);
		WaitingSplashThread th = new WaitingSplashThread () {
			public void run () {
				int count = 0;
				for (Overlap o: overlaps) {
					if (mustStop) return;
					o.autoCuration ();
					if (splash != null) splash.changeStatus (count + " / " + overlaps.size (), count);
					count ++;
				}
				if (splash != null) splash.setVisible (false);
			}
		};
		if (CodonConfiguration.ui) {
			splash = new WaitingSplash (null, "Removing overlaps", "0 / " + overlaps.size () , overlaps.size (), th, true);
			th.start ();
			splash.setVisible (true);
		}
		else {
			th.run ();
		}
		for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) if (i != j) removeSuperpositions (orfsPorFrame [i], orfsPorFrame [j]);
		for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) if (i != j) removeSuperpositions (orfsPorFrame [i], orfsPorFrame [j]);
		for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) if (i != j) removeSuperpositions (orfsPorFrame [i], orfsPorFrame [j]);
		for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) if (i != j) removeSuperpositions (orfsPorFrame [i], orfsPorFrame [j]);
	}
	
	public void removeLowAccuracy (int threehold) {
		for (List <ORF> orfs: orfsPorFrame) {
			for (ORF orf: orfs) {
				if (!orf.isRemoved () && orf.isBlasted ()) {
					if (orf.getAccuracy () < threehold || orf.getSpecificity () < threehold) {
						orf.setRemoved (true);
					}
				}
			}
		}
	}
	
	public StringBuffer wholeSequence = new StringBuffer ("");
	public String amidoSequence [] = new String [6];
	
	public ArrayList <GeneFeature> geneFeatures = new ArrayList <> ();
	
	public void loadBaseSequence () {
		wholeSequence = new StringBuffer ("");
		try {
			BufferedReader reader;
			reader = new BufferedReader (new FileReader (dataSrcFileS));
			String l = reader.readLine ();
			while ((l = reader.readLine ()) != null) {
				wholeSequence.append (l.toUpperCase ());
			}
			reader.close ();
		} catch (FileNotFoundException e) {
			e.printStackTrace ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
		loadAminoSequences ();
	}

	public void loadAminoSequences () {
		String ref = wholeSequence.toString ();
		if (wholeSequence.length () < 3) return;
		amidoSequence [0] = translateFromBaseToAmino (ref, false);
		amidoSequence [3] = translateFromBaseToAmino (ref, true);
		ref = ref.substring (1);
		amidoSequence [1] = translateFromBaseToAmino (ref, false);
		amidoSequence [4] = translateFromBaseToAmino (ref, true);
		ref = ref.substring (1);
		amidoSequence [2] = translateFromBaseToAmino (ref, false);
		amidoSequence [5] = translateFromBaseToAmino (ref, true);
	}
	
	public Vector <Float> getMediaGC (int start, int end) {
		LinkedList <Character> windowGC = new LinkedList <> ();
		Vector <Float> mediaMovelGC = new Vector <> ();
		int sumGC = 0;
		for (int i = start - GC_WINDOW_SIZE; i < Math.min (end, wholeSequence.length ()); i++) {
			if (i < 0) continue;
			char c = wholeSequence.charAt (i);
			windowGC.addLast (c);
			if (c == 'G' || c == 'C') sumGC ++;
			if (windowGC.size () > GC_WINDOW_SIZE) {
				c = windowGC.removeFirst ();
				if (c == 'G' || c == 'C') sumGC --;
			}
			if (i >= start) mediaMovelGC.add ( ((float) sumGC) / windowGC.size ());
		}
		return mediaMovelGC;
	}
	
	static public String translateFromBaseToAmino (String bases, boolean reverse) {
		StringBuffer translation = new StringBuffer ();
		if (reverse) bases = reverseComplementar (bases);
		for (int i = 0; i < bases.length () - 4; i += 3) {
			String codon = bases.substring (i, i + 3);
			translation.append (baseAcidoMap.get (codon));
		}
		return translation.toString ();
	}
       
	public void computeDuplicatedEntries () {
		for (List <ORF> orfs: orfsPorFrame) {
			for (ORF orf: orfs) {
				if (!orf.isRemoved () && orf.isBlasted ()) {
					if (orf.product.blasted && orf.product.hasEntries) duplicatedEntries.addString (orf.product.entryId, orf);
				}
			}
		}
	}
	
	public void computeDuplicatedGenes () {
		for (List <ORF> orfs: orfsPorFrame) {
			for (ORF orf: orfs) {
				if (!orf.isRemoved () && orf.isBlasted ()) {
					if (orf.product.blasted && orf.product.hasEntries && !orf.product.productId.equals ("?")) duplicatedGenes.addString (orf.product.productId, orf);
				}
			}
		}
	}
	
	public static String reverseComplementar (String bases) {
		StringBuffer strBuff = new StringBuffer ();
		for (int i = bases.length () - 1; i >= 0; i --) {
			switch (bases.charAt (i)) {
			case 'A':
				strBuff.append ('T');
				break;
			case 'T':
				strBuff.append ('A');
				break;
			case 'C':
				strBuff.append ('G');
				break;
			case 'G':
				strBuff.append ('C');
				break;
			default:
				break;
			}
		}
		return strBuff.toString ();
	}
	
	public static String complementar (String bases) {
		StringBuffer strBuff = new StringBuffer ();
		for (int i = 0; i < bases.length (); i ++) {
			switch (bases.charAt (i)) {
			case 'A':
				strBuff.append ('T');
				break;
			case 'T':
				strBuff.append ('A');
				break;
			case 'C':
				strBuff.append ('G');
				break;
			case 'G':
				strBuff.append ('C');
				break;
			default:
				break;
			}
		}
		return strBuff.toString ();
	}
	
	public static String reverse (String bases) {
		StringBuffer strBuff = new StringBuffer ();
		for (int i = bases.length () - 1; i >= 0; i --) {
			strBuff.append (bases.charAt (i));
		}
		return strBuff.toString ();
	}
	
	public static String translateFromAminoToBase (String aminos) {
		StringBuffer translation = new StringBuffer ();
		for (int i = 0; i < aminos.length (); i ++) {
			translation.append (acidoBaseMap.get (aminos.charAt (i) + ""));
		}
		return translation.toString ();
	}
	
	public void loadPreannotated () {
		boolean full = true;
		File f = new File (CodonConfiguration.workspace + dataDestDirS + "/fullannotation.txt");
		if (!f.exists ()) {
			f = new File (CodonConfiguration.workspace + dataDestDirS + "/preannotation.txt");
			full = false;
		}
		if (f.exists ()) {
			try {
				ArrayList <ORF> orfs = new ArrayList <ORF> ();						
				BufferedReader reader = new BufferedReader (new FileReader (f));
				String s = null;
				while ((s = reader.readLine ()) != null) {
					String orfS [] = s.split ("\\|");
					ORF orf = new ORF (this, Integer.parseInt (orfS [1]), Integer.parseInt (orfS [2]), Integer.parseInt (orfS [3]));
					orf.id = orfS [0];
					orf.setForced (true);
					orfs.add (orf);
				}
				reader.close ();
				insertOrfs (orfs, true);
				for (ORF orf: orfs) orf.setForced (false);
			} catch (IOException e) {
				e.printStackTrace ();
			}
		}
		if (full) {
			for (ArrayList <ORF> orfs: orfsPorFrame) {
				for (ORF orf: orfs) {
					if (!orf.isBlasted ()) orf.setRemoved (true);
				}
			}
		}
	}
	
	public void annotateRNAOrf (boolean rebuild) {
		File f = new File (CodonConfiguration.workspace + dataDestDirS + "/rna.txt");
		if (!rebuild && f.exists ()) {
			try {
				BufferedReader reader = new BufferedReader (new FileReader (f));
				String s = null;
				while ((s = reader.readLine ()) != null) {
					String orfS [] = s.split ("\\|");
					ORF orf = new ORF (this, Integer.parseInt (orfS [1]), Integer.parseInt (orfS [2]), Integer.parseInt (orfS [3]));
					orf.id = orfS [0];
					orf.setForced (true);
					rnaOrfs.add (orf);
				}
				reader.close ();
			} catch (IOException e) {
				e.printStackTrace ();
			}
		}
		else {
			WaitingSplashThread th = new WaitingSplashThread () {
				public void run () {
					FSMOrfRetriever orfRetriever = new FSMOrfRetriever ();
					ArrayList <ORF> rnas = orfRetriever.retreiveCodonORF (OrfRecognizer.this, "data/bd/tmRNA.fasta", "data/bd/", "tmRNA", 15, 200000, 0, false, true);
					rnaOrfs = rnas;
					if (splash != null) splash.changeStatus ("mRNA", 1);
					rnas = orfRetriever.retreiveCodonORF (OrfRecognizer.this, "data/bd/mRNA.fasta", "data/bd/", "mRNA", 15, 200000, 0, false, true);
					rnaOrfs.addAll (rnas);
					if (splash != null) splash.changeStatus ("miscRNA", 2);
					rnas = orfRetriever.retreiveCodonORF (OrfRecognizer.this, "data/bd/miscRNA.fasta", "data/bd/", "miscRNA", 15, 200000, 0, false, true);
					rnaOrfs.addAll (rnas);
					if (splash != null) splash.changeStatus ("ncRNA", 3);
					rnas = orfRetriever.retreiveCodonORF (OrfRecognizer.this, "data/bd/ncRNA.fasta", "data/bd/", "ncRNA", 15, 200000, 0, false, true);
					rnaOrfs.addAll (rnas);
					if (splash != null) splash.changeStatus ("precursorRNA", 4);
					rnas = orfRetriever.retreiveCodonORF (OrfRecognizer.this, "data/bd/precursorRNA.fasta", "data/bd/", "precursorRNA", 15, 200000, 0, false, true);
					rnaOrfs.addAll (rnas);
					if (splash != null) splash.changeStatus ("rRNA", 5);
					rnas = orfRetriever.retreiveCodonORF (OrfRecognizer.this, "data/bd/rRNA.fasta", "data/bd/", "RNA", 15, 200000, 0, false, true);
					rnaOrfs.addAll (rnas);
					if (splash != null) splash.changeStatus ("tRNA", 6);
					rnas = orfRetriever.retreiveCodonORF (OrfRecognizer.this, "data/bd/tRNA.fasta", "data/bd/", "tRNA", 15, 200000, 0, false, true);
					rnaOrfs.addAll (rnas);
					for (ORF orf: rnaOrfs) {
						if (orf.isForward ()) orf.setOriginalStop (orf.stop);
						else orf.setOriginalStart (orf.start);
					}
					rnaOrfs.sort (new Comparator <ORF> () {
						public int compare (ORF o1, ORF o2) {
							if (o1.start < o2.start) return -1;
							if (o1.start > o2.start) return 1;
							return 0;
						}
					});
					try {
						BufferedWriter writer = new BufferedWriter (new FileWriter (CodonConfiguration.workspace + dataDestDirS + "/rna.txt"));
						for (ORF orf: rnaOrfs) {
							ProductProtein.save (orf.product, CodonConfiguration.workspace + dataDestDirS + "/orfs/" + orf.id + ".unip");
							writer.write (orf.id + "|" + orf.start + "|" + orf.stop + "|" + orf.frame + "|" + orf.isReverse ());
							writer.newLine ();
						}
						writer.close ();
					} catch (IOException e) {
						e.printStackTrace ();
					}
					if (splash != null) splash.setVisible (false);
				}
			};
			if (CodonConfiguration.ui) {
				splash = new WaitingSplash (null, "RNA annotation", "tmRNA", 7, th, true);
				th.start ();
				splash.setVisible (true);
			}
			else {
				th.run ();
			}
		}
		for (ORF orf: rnaOrfs) {
			orf.loadProduct ();
		}
		insertOrfs (rnaOrfs, true);
		for (ORF orf: rnaOrfs) orf.setForced (false);
	}
	
	public void insertOrf (ORF toI, boolean remove) {
		int frameO = toI.frame;
		for (int frame = 0; frame < orfsPorFrame.length; frame++) {
			ArrayList <ORF> orfs = orfsPorFrame [frame];
			for (int i = 0; i < orfs.size (); i ++) {
				ORF orf = orfs.get (i);
				if (orf.stop < toI.start || orf.isForced ()) continue;
				if (orf.stop < toI.stop || orf.start < toI.stop) {
					if (remove) {
						orfs.remove (i);
						i --;
					}
				}
				else {
					if (frame == frameO) {
						orfs.add (i, toI);
					}
					break;
				}
			}
		}
	}
	
	public void insertOrfs (List <ORF> orfsToInsert, boolean remove) {
		for (ORF toI: orfsToInsert) {
			insertOrf (toI, remove);
		}
	}
	
	public void saveAsCodonProject (File file) {
//		int nb = 0;	
//		for (int i = 0; i < orfsPorFrame.length; i++) nb += orfsPorFrame [i].size ();
//		WaitingSplashThread th = new WaitingSplashThread () {
//			public void run () {
//				int nb = 0;	
//				for (int i = 0; i < orfsPorFrame.length; i++) nb += orfsPorFrame [i].size ();
//				int count = 0;
//				for (int i = 0; i < orfsPorFrame.length; i++) {
//					List <ORF> orfs = orfsPorFrame [i];
//					for (ORF orf: orfs) {
//						count ++;
//						if (orf.start != orf.startOriginal || orf.stop != orf.stopOriginal) orf.saveAlterations ();
//						if (splash != null) splash.changeStatus (count + " / " + nb, count);
//					}
//				}
//				if (splash != null) splash.setVisible (false);
//			}
//		};
//		if (CodonConfiguration.ui) {
//			splash = new WaitingSplash (null, "Saving altered ORFs", "0 / " + nb, nb, th, true);
//			th.start ();
//			splash.setVisible (true);
//		}
//		else {
//			th.run ();
//		}
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter (file));
			writer.write ("CODON-V2");
			writer.newLine ();
			writer.write (this.dataDestDirS);
			writer.newLine ();
			writer.write (this.numberA + "|" + this.numberT + "|" + this.numberG + "|" + this.numberC);
			writer.newLine ();
			for (int i = 0; i < wholeSequence.length (); i += 256) {
				if (i + 256 < wholeSequence.length ()) writer.write (wholeSequence.substring (i, i + 256));
				else writer.write (wholeSequence.substring (i));
				writer.newLine ();
			}
			writer.write ("#");writer.newLine ();
			for (int i = 0; i < startStopPorFrame.length; i++) {
				List <Codon> codons = startStopPorFrame [i];
				for (Codon codon: codons) {
					writer.write (codon.toString ());
					writer.newLine ();
				}
				writer.write ("#");writer.newLine ();
			}
			for (int i = 0; i < orfsPorFrame.length; i++) {
				List <ORF> orfs = orfsPorFrame [i];
				for (ORF orf: orfs) {
					if (orf.product != null && orf.product.entryId != null && !orf.product.entryId.equals ("RNA") && !orf.product.entryId.equals ("tRNA")) {
						writer.write (orf.toString ());
						writer.newLine ();
					}
				}
				writer.write ("#");writer.newLine ();
			}
			writer.write ("#");writer.newLine ();
			for (ORF orf: rnaOrfs) {
				writer.write (orf.toString ());
				writer.newLine ();
			}
			writer.write ("#");writer.newLine ();
			writer.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	public boolean loadedFromProject = false;
	public void loadFromCodonProject (File file) {
		loadedFromProject = true;
		clear ();
		try {
			BufferedReader reader = new BufferedReader (new FileReader (file));
			String line = reader.readLine ();
			int version = 1;
			if (line.indexOf ("CODON") == 0) {
				version = Integer.parseInt (line.substring (line.indexOf ("V") + 1));
				dataDestDirS = reader.readLine ();
			}
			else {
				dataDestDirS = line;
			}
			String nbs [] = reader.readLine ().split ("\\|");
			numberA = Integer.parseInt (nbs [0]);
			numberT = Integer.parseInt (nbs [1]);
			numberG = Integer.parseInt (nbs [2]);
			numberC = Integer.parseInt (nbs [3]);
			
			wholeSequence = new StringBuffer ();
			while ((line = reader.readLine ()) != null) {
				if (line.equals ("#")) break;
				wholeSequence.append (line);
			}
			for (int frame = 0; frame < 6; frame++) {
				while ((line = reader.readLine ()) != null) {
					if (line.equals ("#")) break;
					startStopPorFrame [frame].add (Codon.createFromString (line));
				}
			}
			for (int frame = 0; frame < 6; frame++) {
				while ((line = reader.readLine ()) != null) {
					if (line.equals ("#")) break;
					orfsPorFrame [frame].add (ORF.createFromString (this, frame, line, version));
				}
			}
			while ((line = reader.readLine ()) != null) {
				if (line.equals ("#")) break;
			}
			while ((line = reader.readLine ()) != null) {
				if (line.equals ("#")) break;
				rnaOrfs.add (ORF.createFromString (this, -1, line, version));
			}
			insertOrfs (rnaOrfs, false);
			reader.close ();			
		} catch (IOException e) {
			e.printStackTrace ();
		}
		computeDuplicatedEntries ();
		computeDuplicatedGenes ();
		loadAminoSequences ();
	}
	
	public int length () {
		return wholeSequence.length ();
	}
	
	public void exportAsEmbl (File file) {
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter (file));
			int currentIndex [] = new int [] {0, 0, 0, 0, 0, 0};
			List <ORF> orfsPorFrames [] = orfsPorFrame;
			ORF currentOrfs [] = new ORF [] {
					orfsPorFrames [0].size () == 0 ? null : orfsPorFrames [0].get (0), 
					orfsPorFrames [1].size () == 0 ? null : orfsPorFrames [1].get (0), 
					orfsPorFrames [2].size () == 0 ? null : orfsPorFrames [2].get (0), 
					orfsPorFrames [3].size () == 0 ? null : orfsPorFrames [3].get (0), 
					orfsPorFrames [4].size () == 0 ? null : orfsPorFrames [4].get (0), 
					orfsPorFrames [5].size () == 0 ? null : orfsPorFrames [5].get (0)};
			int total = 0;
			for (int i = 0; i < 6; i++) total += orfsPorFrames [i].size ();
			while (total > 0) {
				int min = -1;
				int stMin = Integer.MAX_VALUE;
				for (int i = 0; i < currentOrfs.length; i++) {
					if (currentIndex [i] < orfsPorFrames [i].size () && currentOrfs [i].start < stMin) {
						stMin = currentOrfs [i].start;
						min = i;
					}
				}
				ORF orf = currentOrfs [min];
				if (!orf.isRemoved () && orf.product != null && orf.product.entryId != null && !orf.product.entryId.equals ("RNA") && !orf.product.entryId.equals ("tRNA")) {
					if (min < 3) writer.write ("FT   CDS             " + (orf.start + 1) + ".." + (orf.stop + 3));
					else writer.write ("FT   CDS             complement(" + (orf.start + 1) + ".." + (orf.stop + 3) + ")");
					
					writer.newLine ();
					String s = "/translation=\"" + orf.getAminoSequence () + "\"";
					while (s.length () > 59) {
						writer.write ("FT                   " + s.substring (0, 59));
						s = s.substring (59);
						writer.newLine ();
					}
					if (s.length () > 0) {
						writer.write ("FT                   " + s);
						writer.newLine ();
					}
					orf.loadProduct ();
//					boolean fragment = orf.product.getSpecificity () < 80 && orf.product.identity > 90;
//					String productName = orf.product.identity > 80 ? orf.product.productName : "Uncharacterized Protein"; 
//					if (fragment) productName += " (Fragment)";
					String productName = orf.product.isGene () ? orf.product.gene : orf.product.productName; 
					s = "/product=\"" + productName + "\"";
					while (s.length () > 59) {
						writer.write ("FT                   " + s.substring (0, 59));
						s = s.substring (59);
						writer.newLine ();
					}
					if (s.length () > 0) {
						writer.write ("FT                   " + s);
						writer.newLine ();
					}
					if (orf.isGene () /*&& orf.product.identity > 80*/) {
						writer.write ("FT                   /gene=\"" + orf.product.gene + "\"");
						writer.newLine ();					
					}
					if (orf.product.goTerm.length () > 0) {
						String goT [] = orf.product.goTerm.split ("\\|");
						for (int i = 0; i < goT.length; i++) {
							writer.write ("FT                   /go=\"" + goT [i] + "\"");
							writer.newLine ();
						}
					}
					if (orf.product.kegg.length () > 0) {
						String kegg [] = orf.product.kegg.split ("\\|");
						for (int i = 0; i < kegg.length; i++) {
							writer.write ("FT                   /EC_number=\"" + kegg [i] + "\"");
							writer.newLine ();
						}
					}
				}
				currentIndex [min] ++;
				if (currentIndex [min] < orfsPorFrames [min].size ()) currentOrfs [min] =  orfsPorFrames [min].get (currentIndex [min]);
				total --;
			}
			
			ArrayList <ORF> rnas = new ArrayList <> ();
//			rnas.addAll (tRNAOrfs);
			rnas.addAll (rnaOrfs);
			rnas.sort (new Comparator <ORF> () {
				public int compare (ORF orf0, ORF orf1) {
					if (orf0.start < orf1.start) return -1;
					if (orf0.start > orf1.start) return 1;
					return 0;
				};
			});
			for (ORF orf: rnas) {
				int frame = orf.frame;
				if (orf.product == null || orf.product.entryId == null) continue;
				if (orf.product.entryId.equals ("RNA")) {
					if (frame < 3) writer.write ("FT   RNA             " + (orf.start + 1) + ".." + (orf.stop));
					else writer.write ("FT   RNA             complement(" + (orf.start + 1) + ".." + (orf.stop) + ")");
				}
				else if (orf.product.entryId.equals ("tRNA")) {
					if (frame < 3) writer.write ("FT   tRNA            " + (orf.start + 1) + ".." + (orf.stop));
					else writer.write ("FT   tRNA            complement(" + (orf.start + 1) + ".." + (orf.stop) + ")");
				}
				writer.newLine ();
				String s = "/product=\"" + orf.product.productName + "\"";
				while (s.length () > 59) {
					writer.write ("FT                   " + s.substring (0, 59));
					s = s.substring (59);
					writer.newLine ();
				}
				if (s.length () > 0) {
					writer.write ("FT                   " + s);
					writer.newLine ();
				}
				s = "/translation=\"" + orf.product.matchingSequence + "\"";
				while (s.length () > 59) {
					writer.write ("FT                   " + s.substring (0, 59));
					s = s.substring (59);
					writer.newLine ();
				}
				if (s.length () > 0) {
					writer.write ("FT                   " + s);
					writer.newLine ();
				}
			}
			writer.write ("XX");
			writer.newLine ();
			writer.write ("SQ   Sequence " + length () + " BP; " + numberA + " A; " + numberC + " C; " + numberG + " G; " + numberT + " T; 0 other;");
			writer.newLine ();
			
			int i = 0;
			for (; i < wholeSequence.length (); i += 10) {
				if (i % 60 == 0) {
					if (i != 0) {
						writer.write (" " + i);
						writer.newLine ();
					}
					if (i != wholeSequence.length () - 1) writer.write ("     ");
				}
				writer.write (wholeSequence.substring (i, Math.min (i + 10, wholeSequence.length ())) + " ");
			}
			if (i != wholeSequence.length () - 1) { 
				int k = i % 60;
				int printed = 5 + 11 * (k / 10) + k % 10;
				for (int j = printed; j < 74; j ++) writer.write (" ");
				writer.write (" " + i);
				writer.newLine ();
			}
			
			/*BufferedReader reader = new BufferedReader (new FileReader (dataSrcFileS));
			reader.readLine ();
			String s;
			StringBuffer txt = new StringBuffer ();
			int index = 0;
			while ((s = reader.readLine ()) != null) {
				txt.append (s);
				while (txt.length () > 60) {
					index += 60;
					writer.write ("     ");
					for (int i = 0; i < 6; i++) writer.write (txt.substring (i * 10, (i + 1) * 10) + " ");
					txt = txt.delete (0, 60);
					String sIndex = "" + index;
					while (sIndex.length () < 9) sIndex = " " + sIndex;
					writer.write (sIndex);
					writer.newLine ();
				}
			}
			if (txt.length () > 0) {
				s = "     ";
				index += txt.length ();
				while (txt.length () > 10) {
					s += txt.substring (0, 10) + " ";
					txt.delete (0, 10);
				}
				s +=  txt;
				while (s.length () != 71) s = s + " ";
				writer.write (s);
				String sIndex = "" + index;
				while (sIndex.length () < 9) sIndex = " " + sIndex;
				writer.write (sIndex);
				writer.newLine ();
			}
			reader.close ();*/
			writer.write ("//");
			writer.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	public void loadFromEmbl (File file) {
		clear ();
		this.dataSrcFileS = file.getAbsolutePath ();
		try {
			BufferedReader reader = new BufferedReader (new FileReader (file));
			this.dataDestDirS = file.getName ();
			this.dataDestDirS = dataDestDirS.substring (0, dataDestDirS.indexOf ('.')); 
			File dir = new File (CodonConfiguration.workspace + dataDestDirS + "/orfs");
			if (!dir.exists ()) dir.mkdirs ();
			String line = "";
			String header = "";
			GeneFeature currentGeneFeature = null;
			ORF currentOrf = null;
			int edit = 0;
			clear ();
			while ((line = reader.readLine ()) != null) {
				boolean reverse = false;
				boolean pseudo = false;
				if (line.indexOf ("gene") == 5) {
					if (line.contains (">")) {
						line = line.replace (">", "");
						pseudo = true;
					}
					if (line.contains ("<")) {
						line = line.replace ("<", "");
						pseudo = true;
					}
					if (line.contains ("join")) {
						line = line.replaceAll ("\\d*,\\d*\\.\\.", "").replace ("join(", "").replaceFirst ("\\)", "");
					}
					if (line.contains ("complement")) {
						line = line.replace ("complement(", "");
						line = line.replace (")", "");
						reverse = true;
					}
					int start = Integer.parseInt (line.substring (21, line.indexOf ('.')));
					String s = line.substring (line.indexOf ('.') + 2);
					if (s.contains (",")) s = s.substring (0, s.indexOf (","));
					int stop = Integer.parseInt (s);
					if (start < stop && currentGeneFeature != null) {
						geneFeatures.add (currentGeneFeature);
					}
					currentGeneFeature = new GeneFeature (start - 1, stop - 3, reverse, pseudo); 
					edit = 1;
				}
				else if (line.indexOf ("CDS") == 5) {
					if (line.contains (">")) {
						line = line.replace (">", "");
						pseudo = true;
					}
					if (line.contains ("<")) {
						line = line.replace ("<", "");
						pseudo = true;
					}
					if (line.contains ("join")) {
						line = line.replaceAll ("\\d*,\\d*\\.\\.", "").replace ("join(", "").replaceFirst ("\\)", "");
					}
					if (line.contains ("complement")) {
						line = line.replace ("complement(", "");
						line = line.replace (")", "");
						reverse = true;
					}
					int start = Integer.parseInt (line.substring (21, line.indexOf ('.')));
					String s = line.substring (line.indexOf ('.') + 2);
					if (s.contains (",")) s = s.substring (0, s.indexOf (","));
					int stop = Integer.parseInt (s);
					start -= 1;
					stop -= 3;
					int frame = (start + 1) % 3 + (reverse ? 3: 0);
					
					currentOrf = new ORF (this, start, stop, frame);
					currentOrf.setRemoved (false);
					if (start < stop) {
						orfsPorFrame [frame].add (currentOrf);
					}
					else {
						System.err.println ("Ignored CDS: start = " + start + ", stop = " + stop);
					}
					edit = 2;
				}
				else if (line.indexOf ("RNA ") == 5 || line.indexOf ("tRNA") == 5) {
					if (line.contains (">")) {
						line = line.replace (">", "");
						pseudo = true;
					}
					if (line.contains ("<")) {
						line = line.replace ("<", "");
						pseudo = true;
					}
					if (line.contains ("join")) {
						line = line.replaceAll ("\\d*,\\d*\\.\\.", "").replace ("join(", "").replaceFirst ("\\)", "");
					}
					if (line.contains ("complement")) {
						line = line.replace ("complement(", "");
						line = line.replace (")", "");
						reverse = true;
					}
					int start = Integer.parseInt (line.substring (21, line.indexOf ('.')));
					int stop = Integer.parseInt (line.substring (line.indexOf ('.') + 2));
					start -= 1;
					stop -= 3;
					int frame = (start - 2 ) % 3 + (reverse ? 3: 0);
					currentOrf = new ORF (this, start, stop, frame);
					currentOrf.setRemoved (false);
					rnaOrfs.add (currentOrf);
					if (line.indexOf ("RNA ") == 5)	edit = 3;
					else edit = 4;
				}
				else if (line.indexOf ("SQ") == 0) {
					edit = 5;
				}
				else if (line.indexOf ("ORIGIN") == 0) {
					edit = 5;
				}
				else {
					switch (edit) {
					case 0:
						header += line + "�";
						break;
					case 1:
						if (line.contains ("/gene=")) {
							String gene = line.substring (line.indexOf ("=") + 2, line.length () - 1);
							currentGeneFeature.gene = gene;
						}
						break;
					case 2:
						break;
					case 3:
						if (line.contains ("/product")) {
							String molecule = line.substring (line.indexOf ("/product=") + 11, line.length () - 1);
							currentOrf.product = new ProductProtein ();
							currentOrf.product.entryId = "RNA";
							currentOrf.product.productName = molecule;
							currentOrf.product.productId = molecule;
							currentOrf.product.identity = 100;
							currentOrf.product.setStartMatchSeq (0);
							currentOrf.product.setEndMatchSeq (currentOrf.length / 3);
							currentOrf.product.setStartQuerySeq (0);
							currentOrf.product.setEndQuerySeq (currentOrf.length / 3);
						}
						break;
					case 4:
						if (line.contains ("/product")) {
							String molecule = line.substring (line.indexOf ("/product=") + 11, line.length () - 1);
							currentOrf.product = new ProductProtein ();
							currentOrf.product.entryId = "tRNA";
							currentOrf.product.productName = molecule;
							currentOrf.product.productId = molecule;
							currentOrf.product.identity = 100;
							currentOrf.product.setStartMatchSeq (0);
							currentOrf.product.setEndMatchSeq (currentOrf.length / 3);
							currentOrf.product.setStartQuerySeq (0);
							currentOrf.product.setEndQuerySeq (currentOrf.length / 3);
						}
						break;
					case 5:
						if (line.contains ("//")) break;
						for (int i = 0; i < line.length (); i++) {
							switch (line.charAt (i)) {
								case 'a':
								case 'A':
									wholeSequence.append ('A');
									break;
								case 'T':
								case 't':
									wholeSequence.append ('T');
									break;
								case 'G':
								case 'g':
									wholeSequence.append ('G');
									break;
								case 'C':
								case 'c':
									wholeSequence.append ('C');
									break;
								default:
									break;
							}
						}
						break;
					default:
						break;
					}
				}
			}
			reader.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
		rnaOrfs.sort (new Comparator<ORF> () {
			public int compare (ORF orf0, ORF orf1) {
				if (orf0.start < orf1.start) return -1;
				if (orf0.start == orf1.start) return 0;
				return 1;
			}
		});
//		for (ORF o: tRNAOrfs) {
//			o.product.matchingSequence = translateFromBaseToAmino (o.getSequence (), o.isReverse ());
//			o.product.save (dataDestDirS, o.id, false);
//		}
		for (ORF o: rnaOrfs) {
			o.product.matchingSequence = translateFromBaseToAmino (o.getSequence (), o.isReverse ());
			o.product.save (dataDestDirS, o.id, false);
		}
		insertOrfs (rnaOrfs, false);
//		insertOrfs (tRNAOrfs, false);
		
		OrfRecognizer rnaReco = new OrfRecognizer ();
		for (int i = 0; i < wholeSequence.length (); i++) {
			rnaReco.addChar (wholeSequence.charAt (i));
		}
		startStopPorFrame = rnaReco.startStopPorFrame;
		loadAminoSequences ();
	}
		
	public ArrayList <OrganismFamilly> getOrganisms () {
		ArrayList <OrganismFamilly> famillies = new ArrayList <OrganismFamilly> ();
		for (int i = 0; i < 6; i ++) {
			List <ORF> orfs = orfsPorFrame [i];
			for (ORF orf: orfs) {
				if (!orf.isRemoved () && orf.product != null && orf.product.entryId != null && orf.isBlasted () && !orf.product.entryId.equals ("RNA") && !orf.product.entryId.equals ("tRNA")) {
					boolean found = false;
					for (OrganismFamilly fam: famillies) {
						if (fam.isFromSameFamilly (orf.product.organism)) {
							fam.add (orf);
							found = true;
							break;
						}
					}
					if (!found) {
						OrganismFamilly fam = new OrganismFamilly (orf.product.organism);
						fam.add (orf);
						famillies.add (fam);
					}
				}
			}
		}
		famillies.sort (new Comparator <OrganismFamilly> () {
			public int compare (OrganismFamilly fam0, OrganismFamilly fam1) {
				if (fam0.occurrences > fam1.occurrences) return -1;
				if (fam0.occurrences < fam1.occurrences) return 1;
				return 0;
			}
		});
		for (OrganismFamilly fam: famillies) {
			fam.sort (new Comparator <Organism> () {
				public int compare (Organism o0, Organism o1) {
					if (o0.occurrences > o1.occurrences) return -1;
					if (o0.occurrences < o1.occurrences) return 1;
					return 0;
				};
			});
		}
		return famillies;
	}
	
	public void replace (int selBegin, int selEnd, String newS) {
		wholeSequence = new StringBuffer (wholeSequence.substring (0, selBegin) + newS + wholeSequence.substring (selEnd));
		int diff = selEnd - selBegin - newS.length ();
		
		@SuppressWarnings("unchecked")
		ArrayList <ORF> newOrfsPorFrame [] = new ArrayList [6];
		@SuppressWarnings("unchecked")
		ArrayList <Codon> newStartStopPorFrame [] = new ArrayList [6];  
		
		for (int frame = 0; frame < 6; frame++) {
			newStartStopPorFrame [frame] = new ArrayList <Codon> ();
			ArrayList <Codon> ss = startStopPorFrame [frame];
			ArrayList <Codon> newSs = newStartStopPorFrame [frame];
			for (Codon c: ss) {
				if (c.index < selBegin) {
					newSs.add (c);
				}
			}
		}
		for (int i = Math.max (0, selBegin - 2); i < selBegin + newS.length (); i ++) {
			String codon = wholeSequence.substring (i, i + 3);
			String codonR = reverseComplementar (codon);
			String amino = baseAcidoMap.get (codon);
			String aminoR = baseAcidoMap.get (codonR);
			int frame = (3 + i - 2) % 3;
			int frameR = (3 + i - 2) % 3 + 3;
			if (amino.equals ("M")) newStartStopPorFrame [frame].add (new Codon (START, i));
			if (amino.equals ("V")) newStartStopPorFrame [frame].add (new Codon (VALINA, i));
			if (amino.equals ("L")) newStartStopPorFrame [frame].add (new Codon (LEUCINA, i));
			if (amino.equals ("*")) newStartStopPorFrame [frame].add (new Codon (STOP, i));
			if (amino.equals ("#")) newStartStopPorFrame [frame].add (new Codon (STOP, i));
			if (amino.equals ("+")) newStartStopPorFrame [frame].add (new Codon (STOP, i));
			if (aminoR.equals ("M")) newStartStopPorFrame [frameR].add (new Codon (START_R, i));
			if (aminoR.equals ("V")) newStartStopPorFrame [frameR].add (new Codon (VALINA_R, i));
			if (aminoR.equals ("L")) newStartStopPorFrame [frameR].add (new Codon (LEUCINA_R, i));
			if (aminoR.equals ("*")) newStartStopPorFrame [frameR].add (new Codon (STOP_R, i));
			if (aminoR.equals ("#")) newStartStopPorFrame [frameR].add (new Codon (STOP_R, i));
			if (aminoR.equals ("+")) newStartStopPorFrame [frameR].add (new Codon (STOP_R, i));
		}
		for (int frame = 0; frame < 6; frame++) {
			ArrayList <Codon> ss = startStopPorFrame [frame];
			for (Codon c: ss) {
				if (c.index > selEnd) {
					c.index -= diff;
					int newFrame = (c.index - 2) % 3 + 3 * (frame > 2 ? 1 : 0); 
					newStartStopPorFrame [newFrame].add (c);
				}
			}
		}
		
		for (int frame = 0; frame < 6; frame++) {
			newOrfsPorFrame [frame] = new ArrayList <ORF> ();
			ArrayList <ORF> orfs = orfsPorFrame [frame];
			ArrayList <ORF> newOrfs = newOrfsPorFrame [frame];
			for (ORF orf: orfs) {
				if (orf.start < selBegin && orf.stop < selBegin) {
					newOrfs.add (orf);
				}
			}
		}
		for (int frame = 0; frame < 6; frame++) {
			ArrayList <ORF> orfs = orfsPorFrame [frame];
			for (ORF orf: orfs) {
				if (orf.start > selEnd) {
					orf.start -= diff;
					orf.startOriginal -= diff;
					orf.stop -= diff;
					orf.stopOriginal -= diff;
					int newFrame = (orf.start - 2) % 3 + 3 * (orf.isReverse () ? 1 : 0);
					orf.frame = newFrame;
					newOrfsPorFrame [newFrame].add (orf);
				}
			}
		}
		
		orfsPorFrame = newOrfsPorFrame;
		startStopPorFrame = newStartStopPorFrame;
		
		for (int j = 0; j < 6; j++) {
			for (int i = Math.max (0, selBegin - 2); i < selBegin + newS.length (); i += 3) {
				ORF orf = createOrf (j, i);
				if (orf != null) {
					i = orf.stop + 1;
					if (j > 2) {
						i = startStopPorFrame [j].get (getCodonIndex (startStopPorFrame [j], i) + 1).index;
					}
				}
			}
		}
		loadAminoSequences ();
	}
	
	public ORF createOrf (int frame, int index) {
		ArrayList <Codon> ss = startStopPorFrame [frame];
		ArrayList <ORF> orfs = orfsPorFrame [frame];
		int ind = getCodonIndex (ss, index);
		Codon c = ss.get (ind);
		if (frame < 3) {
			while (c.type == OrfRecognizer.STOP) {
				ind ++;
				if (ind > ss.size () - 1) return null;
				c = ss.get (ind);
			}
			while (c.type != OrfRecognizer.STOP) {
				ind --;
				if (ind < 0) break;
				c = ss.get (ind);
			}
			ind ++;
			c = ss.get (ind);
			Codon startC = c;
			while (c.type != OrfRecognizer.STOP) {
				ind ++;
				if (ind > ss.size () - 1) return null;
				c = ss.get (ind);
			}
			Codon stopC = c;
			if (stopC.index - startC.index < CodonConfiguration.min_orf_length) return null;
			ORF orfN = new ORF (this, startC.index, stopC.index, false);
			for (ORF orf: orfs) {
				if (orf.intersection (orfN) != 0) {
					if (orf.isRemoved ()) {
						int i = orfs.indexOf (orf);
						orfs.remove (i);
						orfs.add (i, orfN);
						return orfN;
					}
					return orf;
				}
			}
			for (int i = 0; i < orfs.size (); i ++) {
				if (orfs.get (i).start > orfN.start) {
					orfs.add (i, orfN);
					return orfN;
				}
			}
		}
		else {
			while (c.type != OrfRecognizer.STOP_R) {
				ind --;
				if (ind < 0) return null;
				c = ss.get (ind);
			}
			Codon startC = c;
			ind ++;
			if (ind > ss.size () - 1) return null;
			c = ss.get (ind);
			while (c.type != OrfRecognizer.STOP_R) {
				ind ++;
				if (ind > ss.size () - 1) return null;
				c = ss.get (ind);
			}
			ind --;
			if (ind < 0) return null;
			c = ss.get (ind);
			Codon stopC = c;
			if (stopC.index - startC.index < CodonConfiguration.min_orf_length) return null;
			ORF orfN = new ORF (this, startC.index, stopC.index, true);
			for (ORF orf: orfs) {
				if (orf.intersection (orfN) != 0) {
					if (orf.isRemoved ()) {
						int i = orfs.indexOf (orf);
						orfs.remove (i);
						orfs.add (i, orfN);
						return orfN;
					}
					return orf;
				}
			}
			for (int i = 0; i < orfs.size (); i ++) {
				if (orfs.get (i).start > orfN.start) {
					orfs.add (i, orfN);
					return orfN;
				}
			}
		}
		return null;
	}
	public int getOrfCount () {
		int total = 0;
		for (List <ORF> orfs : orfsPorFrame) for (ORF orf: orfs) if (!orf.isRemoved ()) total ++;
		return total;
	}
	
	public void upgradeRNAS () {
		try {
			Process proc = Runtime.getRuntime ().exec ("java "
					+ "-cp \"lib\\activation-1.1.jar;lib\\animal-sniffer-annotations-1.14.jar;lib\\antlr4-4.5.3.jar;lib\\aopalliance-1.0.jar;lib\\aopalliance-repackaged-2.4.0-b10.jar;lib\\avro-1.7.6.jar;lib\\avro-ipc-1.7.6.jar;lib\\axis2-1.6.2.jar;lib\\checker-compat-qual-2.0.0.jar;lib\\classmate-0.8.0.jar;lib\\commons-beanutils-1.9.2.jar;lib\\commons-collections-3.2.1.jar;lib\\commons-compress-1.12.jar;lib\\commons-io-2.4.jar;lib\\commons-lang-2.6.jar;lib\\commons-lang3-3.6.jar;lib\\commons-logging-1.2.jar;lib\\commons-math3-3.6.1.jar;lib\\commons-net-3.6.jar;lib\\commons-text-1.8.jar;lib\\error_prone_annotations-2.1.3.jar;lib\\failsafe-1.0.3.jar;lib\\FastInfoset-1.2.13.jar;lib\\guava-24.1.1-jre.jar;lib\\guice-4.0.jar;lib\\guice-multibindings-4.0.jar;lib\\hibernate-validator-5.0.1.Final.jar;lib\\hk2-api-2.4.0-b10.jar;lib\\hk2-locator-2.4.0-b10.jar;lib\\hk2-utils-2.4.0-b10.jar;lib\\httpclient-4.5.3.jar;lib\\httpcore-4.4.6.jar;lib\\httpmime-4.5.3.jar;lib\\istack-commons-runtime-3.0.5.jar;lib\\j2objc-annotations-1.1.jar;lib\\jackson-annotations-2.9.5.jar;lib\\jackson-core-2.9.5.jar;lib\\jackson-core-asl-1.9.13.jar;lib\\jackson-databind-2.9.4.jar;lib\\jackson-mapper-asl-1.9.13.jar;lib\\japi-1.0.30.jar;lib\\javaparser-1.0.11.jar;lib\\javassist-3.18.1-GA.jar;lib\\javax.annotation-api-1.2.jar;lib\\javax.inject-1.jar;lib\\javax.inject-2.4.0-b10.jar;lib\\javax.mail-1.6.1.jar;lib\\javax.mail-api-1.6.1.jar;lib\\javax.ws.rs-api-2.0.1.jar;lib\\jaxb2-basics-1.11.1.jar;lib\\jaxb2-basics-runtime-1.11.1.jar;lib\\jaxb2-basics-tools-1.11.1.jar;lib\\jaxb-api-2.3.0.jar;lib\\jaxb-core-2.3.0.jar;lib\\jaxb-impl-2.3.0.jar;lib\\jaxb-runtime-2.3.0.jar;lib\\jboss-logging-3.1.1.GA.jar;lib\\jcl-over-slf4j-1.7.7.jar;lib\\jcommander-1.32.jar;lib\\jersey-client-2.17.jar;lib\\jersey-common-2.17.jar;lib\\jersey-guava-2.17.jar;lib\\jersey-media-jaxb-2.17.jar;lib\\jetty-6.1.26.jar;lib\\jetty-util-6.1.26.jar;lib\\jsr305-1.3.9.jar;lib\\log4j-api-2.12.1.jar;lib\\log4j-core-2.12.1.jar;lib\\logback-classic-1.2.3.jar;lib\\logback-core-1.2.3.jar;lib\\metrics-core-3.2.2.jar;lib\\netty-3.4.0.Final.jar;lib\\noggit-0.8.jar;lib\\osgi-resource-locator-1.0.1.jar;lib\\paranamer-2.3.jar;lib\\servlet-api-2.5-20081211.jar;lib\\slf4j-api-1.7.25.jar;lib\\snappy-java-1.0.5.jar;lib\\solr-solrj-7.2.0.jar;lib\\stax2-api-3.1.4.jar;lib\\stax-ex-1.7.8.jar;lib\\txw2-2.3.0.jar;lib\\validation-api-1.1.0.Final.jar;lib\\velocity-1.7.jar;lib\\woodstox-core-asl-4.4.1.jar;lib\\zookeeper-3.4.10.jar;bin\" -Dswing.defaultlaf=javax.swing.plaf.nimbus.NimbusLookAndFeel"
					+ " codon.UpgradeRNA");
			PrintStream ouP = null;
			ouP = new PrintStream (proc.getOutputStream ());
			ouP.println ("##");
			for (ORF orf: rnaOrfs) {
				ouP.println (orf.getSequence ());
				ouP.println (orf.product.productName);
			}
			ouP.println ("##");
			ouP.flush ();

			InputStream inP = proc.getInputStream ();
			int c;
			while ((c = inP.read ()) != -1) {
				System.out.print ((char) c);
			}
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
}
