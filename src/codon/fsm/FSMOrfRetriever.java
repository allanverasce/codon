package codon.fsm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.zip.GZIPInputStream;

import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.data.ProductProtein;
import codon.recognizer.OrfRecognizer;
import splash.WaitingSplash;
import splash.WaitingSplashThread;

// http://gtrnadb2009.ucsc.edu/download.html

public class FSMOrfRetriever {
	
	HashFSMBuilder builder = null;
	
	public static char alphabetNucleotide [] = new char [] {'A', 'T', 'G', 'C'};
	public static char alphabetAmino [] = new char [] {
		'F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y', '#', '+', 
		'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', '*', 'G', 'X', 'U'
	};
	
	public int buildFSMs (String cdsFilename, String fsmDestDir, int prefixMaxSize, int sequencesPerFSM, int cutStartCodon, char [] alphabet) {
		builder = new HashFSMBuilder ();
		
		File cdsFile = new File (cdsFilename);
		String fName = cdsFile.getName ();
		String destFilename = fsmDestDir + fName.substring (0, fName.lastIndexOf (".") + 1);
		File d = new File (destFilename).getParentFile ();
		
		File lstF [] = d.listFiles ();
		System.out.println(d.getAbsolutePath());
		int machineNb = 0;
		for (File f: lstF) {
			if (f.getAbsolutePath ().contains (destFilename + prefixMaxSize + ".fsm.")) machineNb ++;
		}
		if (machineNb != 0) return machineNb;
		long init = System.currentTimeMillis ();
		try {
			BufferedReader reader = null;
			GZIPInputStream in = null;
			if (cdsFilename.indexOf (".fasta.gz") == cdsFilename.length () - 9) {
				in = new GZIPInputStream (new FileInputStream (cdsFilename));
				reader = new BufferedReader (new InputStreamReader (in));
			}
			else {
				reader = new BufferedReader (new FileReader (cdsFilename));
			}
			String molecule;
			String cds = null;
			int count = 0;
			ArrayList <String> cdsSequences = new ArrayList <String> ();
			ArrayList <String> moleculas = new ArrayList <String> ();
			molecule = reader.readLine ();
			String l;
			while (molecule != null) {
				molecule = molecule.substring (1);
				cds = reader.readLine ();
				moleculas.add (molecule);
				while ((l = reader.readLine ()) != null) {
					if (l.charAt (0) == '>') {
						molecule = l;
						break;
					}
					else {
						cds += l;
					}
				}
				if (l == null) molecule = null;
				if (cutStartCodon != 0) cds = cds.substring (cutStartCodon);
				cdsSequences.add (cds);
				count ++;
				if (count == sequencesPerFSM) {
					builder = new HashFSMBuilder (alphabet);
					builder.create (cdsSequences, moleculas, destFilename + prefixMaxSize + ".fsm." + machineNb, destFilename + prefixMaxSize + ".lex." + machineNb, prefixMaxSize);
					machineNb ++;
					cdsSequences.clear ();
					moleculas.clear ();
					count = 0;
				}
			}
			reader.close ();
			if (in != null) in.close ();
			builder = new HashFSMBuilder (alphabet);
			builder.create (cdsSequences, moleculas, destFilename + prefixMaxSize + ".fsm." + machineNb, destFilename + prefixMaxSize + ".lex." + machineNb, prefixMaxSize);
			System.out.println ("Building time (ms):" + (System.currentTimeMillis () - init));
			return machineNb + 1;
		} catch (IOException e) {
			e.printStackTrace ();
		}
		return -1;
	}
	
	public ORF createORFBases (OrfRecognizer reco, String entryID, int start, String cdsSeq, String molecule, boolean reverse, boolean save) {
		ORF orf = new ORF (reco, start, start + cdsSeq.length () - 3, reverse);
		orf.product = new ProductProtein ();
		orf.product.entryId = entryID;
		if (molecule.contains ("/molecule=")) {
			orf.product.productName = molecule.substring (molecule.indexOf ("=") + 1, molecule.indexOf ("/s") - 1);
		}
		else if (molecule.contains (" ")) {
			orf.product.productName = molecule.substring (0, molecule.indexOf (" "));
		}
		else {
			orf.product.productName = molecule;
		}
		orf.product.productId = entryID;
		orf.product.identity = 100;
		orf.product.setStartMatchSeq (0);
		orf.product.setEndMatchSeq (cdsSeq.length () / 3);
		orf.product.setStartQuerySeq (0);
		orf.product.setEndQuerySeq (cdsSeq.length () / 3);
		orf.product.matchingSequence = OrfRecognizer.translateFromBaseToAmino (cdsSeq, reverse);
//		orf.setForced (true);
		if (save) orf.product.save (reco.dataDestDirS, orf.id, false);
		return orf;
	}
	
	public ORF createORFAmino (OrfRecognizer reco, String entryID, int start, String aminoSeq, String molecule, boolean reverse, boolean save) {
		ORF orf = new ORF (reco, start, start + aminoSeq.length () * 3, reverse);
		orf.product = new ProductProtein ();
		orf.product.entryId = entryID;
		orf.product.productId = entryID;
		if (molecule.contains ("/molecule=")) {
			orf.product.productName = molecule.substring (molecule.indexOf ("=") + 1, molecule.indexOf ("/s") - 1);
		}
		else if (molecule.contains ("tr|") || molecule.contains ("p|")) {
			molecule = molecule.substring (molecule.indexOf ("|") + 1);
			orf.product.entryId = molecule.substring (0, molecule.indexOf ("|"));
			molecule = molecule.substring (molecule.indexOf ("|") + 1);
			orf.product.productId = molecule.substring (0, molecule.indexOf (" ") - 1);
			molecule = molecule.substring (molecule.indexOf (" ") + 1);
			if (molecule.contains ("=")) {
				orf.product.productName = molecule.substring (0, molecule.indexOf ("=") - 3);
				molecule = molecule.substring (molecule.indexOf ("=") - 3);
			}
			else {
				orf.product.productName = molecule;
			}
			if (molecule.contains ("OS=")) {
				molecule = molecule.substring (molecule.indexOf ("OS=") + 3);
				if (molecule.contains ("=")) orf.product.organism = molecule.substring (0, molecule.indexOf ("=") - 3);
				else orf.product.organism = molecule;
			}
			if (molecule.contains (" GN=")) {
				molecule = molecule.substring (molecule.indexOf (" GN=") + 4);
				if (molecule.contains (" ")) orf.product.gene = molecule.substring (0, molecule.indexOf (" "));
				else orf.product.gene = molecule;
				if (orf.product.gene.length () > 7) orf.product.gene = "";
			}
		}
		else if (molecule.contains (" ")) {
			orf.product.productName = molecule.substring (1, molecule.indexOf (" "));
		}
		else {
			orf.product.productName = molecule;
		}
		
		orf.product.identity = 100;
		orf.product.setStartMatchSeq (0);
		orf.product.setEndMatchSeq (aminoSeq.length () + 1);
		orf.product.setStartQuerySeq (0);
		orf.product.setEndQuerySeq (aminoSeq.length () + 1);
		orf.product.matchingSequence = "M" + aminoSeq;
//		orf.setForced (true);
		if (save) orf.product.save (reco.dataDestDirS, orf.id, false);
		return orf;
	}
	
	WaitingSplash splash;
	String forward;
	String reverse;
	
	public ArrayList <ORF> mark (OrfRecognizer reco, String entryId, String refs[], boolean removeOverlaps, boolean searchAmino) {
		ArrayList <ORF> cdsOrfs = new ArrayList <ORF> ();
		for (int k = 0; k < refs.length; k++) {
			ArrayList <RecognizedPattern> marks = builder.mark (refs [k], removeOverlaps);
			for (int i = 0; i < marks.size (); i++) {
				RecognizedPattern mark = marks.get (i);
				String seq = refs [k].substring (mark.start, mark.end);
				String req = mark.recognized;
				if (searchAmino) {
					mark.start = mark.start * 3 + k;
					mark.end = mark.end * 3 + k;
					if (k > 2) {
						mark.start = reco.length () - mark.end;
					}
					cdsOrfs.add (createORFAmino (reco, entryId, mark.start, seq, req, k > 2, false));
				}
				else {
					if (k == 1) mark.start = reco.length () - mark.end;
					cdsOrfs.add (createORFBases (reco, entryId, mark.start, seq, req, k == 1, false));	
				}
			}
		}
		return cdsOrfs;
	}

	public void retreiveCodonORFNonPersistent (OrfRecognizer reco, String cdsFilename, String entryId, int prefixMaxSize, int sequencesPerFSM, int cutStartCodon, boolean searchAmino, FSMOrfRetrieverListener listener) {
		builder = new HashFSMBuilder ();

		File cdsFile = new File (cdsFilename);
		String fName = cdsFile.getName ();
		
		forward = reco.wholeSequence.toString ();
		reverse = OrfRecognizer.reverseComplementar (reco.wholeSequence.toString ());
		
		String refs [];
		if (searchAmino) {
			refs = new String [6];
			refs [0] = OrfRecognizer.translateFromBaseToAmino (forward, false); 
			refs [1] = OrfRecognizer.translateFromBaseToAmino (forward.substring (1), false);
			refs [2] = OrfRecognizer.translateFromBaseToAmino (forward.substring (2), false);
			refs [3] = OrfRecognizer.translateFromBaseToAmino (reverse, false);
			refs [4] = OrfRecognizer.translateFromBaseToAmino (reverse.substring (1), false);
			refs [5] = OrfRecognizer.translateFromBaseToAmino (reverse.substring (2), false);
		}
		else {
			refs = new String [] {forward, reverse};
		}
		
		char [] alphabet = searchAmino ? alphabetAmino : alphabetNucleotide;
		int machineNb = 0;
		try {
			BufferedReader reader = null;
			if (cdsFilename.indexOf (".fasta.gz") == cdsFilename.length () - 9) {
				GZIPInputStream in = new GZIPInputStream (new FileInputStream (cdsFilename));
				reader = new BufferedReader (new InputStreamReader (in));
			}
			else {
				reader = new BufferedReader (new FileReader (cdsFilename));
			}
			String molecule;
			String cds = null;
			int count = 0;
			ArrayList <String> cdsSequences = new ArrayList <String> ();
			ArrayList <String> moleculas = new ArrayList <String> ();
			molecule = reader.readLine ();
			String l;
			int total = 0;
			while (molecule != null) {
				molecule = molecule.substring (1);
				cds = reader.readLine ();
				moleculas.add (molecule);
				while ((l = reader.readLine ()) != null) {
					if (l.length() == 0) System.out.println (molecule + " " + cds);
					if (l.charAt (0) == '>') {
						molecule = l;
						break;
					}
					else {
						cds += l;
					}
				}
				if (l == null) molecule = null;
				if (cutStartCodon != 0) cds = cds.substring (cutStartCodon);
				cdsSequences.add (cds);
				count ++;
				if (count == sequencesPerFSM) {
					System.out.println (cdsFile.getCanonicalPath () + " FSM " + machineNb);
					builder = new HashFSMBuilder (alphabet);
					builder.create (cdsSequences, moleculas, null, null, prefixMaxSize);
					ArrayList <ORF> marked = mark (reco, entryId, refs, false, searchAmino);
					listener.found (marked);
					total += marked.size ();
					System.out.println ("Found " + marked.size () + " total " + total);
					machineNb ++;
					cdsSequences.clear ();
					moleculas.clear ();
					count = 0;
				}
			}
			reader.close ();
			System.out.println (fName + " FSM " + machineNb);
			builder = new HashFSMBuilder (alphabet);
			builder.create (cdsSequences, moleculas, null, null, prefixMaxSize);
			ArrayList <ORF> marked = mark (reco, entryId, refs, false, searchAmino);
			listener.found (marked);
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	public ArrayList <ORF> retreiveCodonORFNonPersistent (OrfRecognizer reco, String cdsFilename, String entryId, int prefixMaxSize, int sequencesPerFSM, int cutStartCodon, boolean searchAmino, boolean removeOverlaps) {
		builder = new HashFSMBuilder ();

		File cdsFile = new File (cdsFilename);
		String fName = cdsFile.getName ();
		
		ArrayList <ORF> cdsOrfs = new ArrayList <ORF> ();
		
		forward = reco.wholeSequence.toString ();
		reverse = OrfRecognizer.reverseComplementar (reco.wholeSequence.toString ());
		
		String refs [];
		if (searchAmino) {
			refs = new String [6];
			refs [0] = OrfRecognizer.translateFromBaseToAmino (forward, false); 
			refs [1] = OrfRecognizer.translateFromBaseToAmino (forward.substring (1), false);
			refs [2] = OrfRecognizer.translateFromBaseToAmino (forward.substring (2), false);
			refs [3] = OrfRecognizer.translateFromBaseToAmino (reverse, false);
			refs [4] = OrfRecognizer.translateFromBaseToAmino (reverse.substring (1), false);
			refs [5] = OrfRecognizer.translateFromBaseToAmino (reverse.substring (2), false);
		}
		else {
			refs = new String [] {forward, reverse};
		}
		
		char [] alphabet = searchAmino ? alphabetAmino : alphabetNucleotide;
		int machineNb = 0;
		try {
			BufferedReader reader = null;
			if (cdsFilename.indexOf (".fasta.gz") == cdsFilename.length () - 9) {
				GZIPInputStream in = new GZIPInputStream (new FileInputStream (cdsFilename));
				reader = new BufferedReader (new InputStreamReader (in));
			}
			else {
				reader = new BufferedReader (new FileReader (cdsFilename));
			}
			String molecule;
			String cds = null;
			int count = 0;
			ArrayList <String> cdsSequences = new ArrayList <String> ();
			ArrayList <String> moleculas = new ArrayList <String> ();
			molecule = reader.readLine ();
			String l;
			while (molecule != null) {
				molecule = molecule.substring (1);
				cds = reader.readLine ();
				moleculas.add (molecule);
				while ((l = reader.readLine ()) != null) {
					if (l.length() == 0) System.out.println (molecule + " " + cds);
					if (l.charAt (0) == '>') {
						molecule = l;
						break;
					}
					else {
						cds += l;
					}
				}
				if (l == null) molecule = null;
				if (cutStartCodon != 0) cds = cds.substring (cutStartCodon);
				cdsSequences.add (cds);
				count ++;
				if (count == sequencesPerFSM) {
					System.out.println (cdsFile.getCanonicalPath () + " FSM " + machineNb);
					builder = new HashFSMBuilder (alphabet);
					builder.create (cdsSequences, moleculas, null, null, prefixMaxSize);
					ArrayList <ORF> marked = mark (reco, entryId, refs, removeOverlaps, searchAmino); 
					cdsOrfs.addAll (marked);
					System.out.println ("Found " + marked.size () + " total " + cdsOrfs.size ());
					machineNb ++;
					cdsSequences.clear ();
					moleculas.clear ();
					count = 0;
				}
			}
			reader.close ();
			System.out.println (fName + " FSM " + machineNb);
			builder = new HashFSMBuilder (alphabet);
			builder.create (cdsSequences, moleculas, null, null, prefixMaxSize);
			cdsOrfs.addAll (mark (reco, entryId, refs, removeOverlaps, searchAmino));
		} catch (IOException e) {
			e.printStackTrace ();
		}

		if (removeOverlaps) {
			ArrayList <ORF> finalORFs = removeOverlaps (cdsOrfs);
			return finalORFs;
		}
		else {
			return cdsOrfs;
		}
	}
	
	long loadTime = 0;
	long markTime = 0;
	
	public ArrayList <ORF> removeOverlaps (ArrayList <ORF> cdsOrfs) {
		ArrayList <ORF> finalORFs = new ArrayList <ORF> ();
		for (int i = 0; i < cdsOrfs.size (); i++) {
			ORF p1 = cdsOrfs.get (i);
			boolean addToFinal = true;
			for (int j = cdsOrfs.size () - 1; j > 0 ; j --) {
				ORF p2 = cdsOrfs.get (j);
				if (p1 != p2 && p1.start == p2.start && p1.stop == p2.stop) {
					cdsOrfs.remove (p2);
				}
			}
			for (int j = 0; j < cdsOrfs.size (); j++) {
				ORF p2 = cdsOrfs.get (j);
				if (p1 != p2) {
					if (p1.start >= p2.start && p1.start < p2.stop || p2.start >= p1.start && p2.start < p1.stop) {
						if (p1.length < p2.length) { 
							addToFinal = false;
							break;
						}
					}
				}
			}
			if (addToFinal) {
				finalORFs.add (p1);
			}
		}
		Collections.sort (finalORFs, new Comparator <ORF>() {
			public int compare (ORF o0, ORF o1) {
				if (o0.start < o1.start) return -1;
				if (o0.start > o1.start) return 1;
				if (o0.start == o1.start) {
					if (o0.stop < o1.stop) return -1;
					if (o0.stop > o1.stop) return 1;
					if (o0.stop == o1.stop) return 0;
					return 0;
				}
				return 0;
			}
		});
		return finalORFs;
	}
	
	public ArrayList <ORF> retreiveCodonORF (OrfRecognizer reco, String cdsFilename, String fsmDestDir, String entryId, int prefixMaxSize, int sequencesPerFSM, int cutStartCodon, boolean searchAmino, boolean removeOverlaps) {
		File cdsFile = new File (cdsFilename);
		String fName = cdsFile.getName ();
		String destFilename = fsmDestDir + fName.substring (0, fName.lastIndexOf (".") + 1);
		
		ArrayList <ORF> cdsOrfs = new ArrayList <ORF> ();
		
		int machineNB = 0;
		machineNB = buildFSMs (cdsFilename, fsmDestDir, prefixMaxSize, sequencesPerFSM, cutStartCodon, searchAmino ? alphabetAmino : alphabetNucleotide);
		forward = reco.wholeSequence.toString ();
		reverse = OrfRecognizer.reverseComplementar (reco.wholeSequence.toString ());
		
		String refs [];
		if (searchAmino) {
			refs = new String [6];
			refs [0] = OrfRecognizer.translateFromBaseToAmino (forward, false); 
			refs [1] = OrfRecognizer.translateFromBaseToAmino (forward.substring (1), false);
			refs [2] = OrfRecognizer.translateFromBaseToAmino (forward.substring (2), false);
			refs [3] = OrfRecognizer.translateFromBaseToAmino (reverse, false);
			refs [4] = OrfRecognizer.translateFromBaseToAmino (reverse.substring (1), false);
			refs [5] = OrfRecognizer.translateFromBaseToAmino (reverse.substring (2), false);
		}
		else {
			refs = new String [] {forward, reverse};
		}
		long stepTime = System.currentTimeMillis ();	
		for (int m = 0; m < machineNB; m ++) {
			System.out.println (cdsFilename + ", FSM " + (m + 1) + "/" + (machineNB));
			builder.load (destFilename + prefixMaxSize + ".fsm." + m, destFilename + prefixMaxSize + ".lex." + m);
			cdsOrfs.addAll (mark (reco, entryId, refs, removeOverlaps, searchAmino));
		}
		if (removeOverlaps) {
			ArrayList <ORF> finalORFs = removeOverlaps (cdsOrfs);
			markTime += System.currentTimeMillis () - stepTime;
			return finalORFs;
		}
		else {
			markTime += System.currentTimeMillis () - stepTime;
			return cdsOrfs;
		}
	}
	
}
