package codon.fsm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;

import splash.WaitingSplash;

public class HashFSMBuilder {
	char [] alphabet = new char [] {'A', 'T', 'G', 'C'};
	
	public HashFSMBuilder () {}
	public HashFSMBuilder (char [] alphabet) {
		this.alphabet = alphabet;
	}
	
	class Transition {
		char read;
		FSMState nextState;
		public Transition (char read, FSMState nextState) {
			super ();
			this.read = read;
			this.nextState = nextState;
		}
	}
	
	static long nbstates = 0;
	
	class FSMState {
		ArrayList <Transition> transistions = new ArrayList <> ();
		String recognized = null;
		long id;
		int count = 1;
		public FSMState () {
			id = nbstates;
			nbstates ++;
		}
		FSMState addTransition (char c) {
			for (Transition tr: transistions) {
				if (tr.read == c) {
					tr.nextState.count ++;
					return tr.nextState;
					
				}
			}
			Transition tr = new Transition (c, new FSMState ());
			transistions.add (tr);
			return tr.nextState;
		}
		FSMState next (char c) {
			for (Transition tr: transistions) {
				if (tr.read == c) return tr.nextState;
			}
			return null;
		}
		void close (String prefix) {
			for (int i = 0; i < alphabet.length; i++) {
				boolean fnd = false;
				for (Transition tr: transistions) {
					if (tr.read == alphabet [i]) { 
						fnd = true;
						break;
					}
				}
				if (!fnd) {
					String pref = prefix + alphabet [i];
					while (pref.length () != 0) {
						pref = pref.substring (1, pref.length ());
						FSMState st = prefixes.get (pref);
						if (st != null) {
							transistions.add (new Transition (alphabet [i], st));
							fnd = true;
							break;
						}
					}
					if (!fnd) transistions.add (new Transition (alphabet [i], idleState));
				}
			}
		}
	}
	
	FSMState idleState = new FSMState ();
	Hashtable <String, FSMState> prefixes = new Hashtable <String, HashFSMBuilder.FSMState> ();
	WaitingSplash splash;

	Hashtable <String, ArrayList <String>> molecules = new Hashtable <> ();
	Hashtable <String, ArrayList <String>> lexique = new Hashtable <> ();
	
	public void create (List <String> strs, List <String> recognized, String fsmFilename, String lexqueFileName, int maxPrefixSize) {
		nbstates = 0;
		idleState = new FSMState ();
		prefixes.clear ();
		prefixes.put ("", idleState);
		
		lexique.clear ();
		molecules.clear ();
		
		int maxSize = 0;
		Hashtable <String, FSMState> workingHash = new Hashtable <> ();
		Hashtable <String, FSMState> finalHash = new Hashtable <> ();
		
		int k = 0;
		for (String s: strs) {
			if (s.length () > maxSize) maxSize=s.length ();
			workingHash.put (s, idleState);
			ArrayList <String> mols = molecules.get (s);
			if (mols == null) {
				mols = new ArrayList <String> ();
				molecules.put (s, mols);
			}
			mols.add (recognized.get (k));
			k ++;
		}
		
		for (int i = 0; i < maxSize; i ++) {
			Enumeration <String> keys = workingHash.keys ();
			while (keys.hasMoreElements ()) {
				String key = keys.nextElement ();
				char transition = key.charAt (i);
				FSMState state = workingHash.get (key);
				state = state.addTransition (transition);
				if (state.count == 1) {
					prefixes.put (key.substring (0, i + 1), state);
				}
				if (i == key.length () - 1) {
					if (i > 30) {
						finalHash.put (key, state);
						state.recognized = key;
					}
					workingHash.remove (key);
				}
				else {
					workingHash.put (key, state);
				}
			}
			if (i > 30) {
				keys = workingHash.keys ();
				while (keys.hasMoreElements ()) {
					String key = keys.nextElement ();
					FSMState state = workingHash.get (key);
					if (state.count == 1) {
						finalHash.put (key, state);
						workingHash.remove (key);
						state.recognized = key.substring (0, i + 1);
					}
				}
			}
			if (i == maxPrefixSize) {
				keys = workingHash.keys ();
				while (keys.hasMoreElements ()) {
					String key = keys.nextElement ();
					FSMState state = workingHash.get (key);
					finalHash.put (key, state);
					state.recognized = key.substring (0, i + 1);
				}
				workingHash.clear ();
			}
			if (workingHash.size () == 0) break;
		}
		
		try {
			BufferedWriter writer = null;
			if (lexqueFileName != null) writer = new BufferedWriter (new FileWriter (lexqueFileName));
			Enumeration <String> keys = finalHash.keys ();
			while (keys.hasMoreElements ()) {
				String key = keys.nextElement ();
				FSMState current = idleState;
				FSMState next;
				char c = key.charAt (0);
				int pos = 0;
				while ((next = current.next (c)) != null) {
					current = next;
					pos ++;
					if (pos == key.length ()) break;
					c = key.charAt (pos);
				}
				String pre = key.substring (0, pos);
				String suf = key.substring (pos);
				ArrayList <String> sufs = lexique.get (pre);
				if (sufs == null) {
					sufs = new ArrayList <String> ();
					lexique.put (pre, sufs);
				}
				sufs.add (suf);				
				if (lexqueFileName != null) {
					ArrayList <String> mols = molecules.get (key);
					String m = "";
					for (String mol: mols) m += ";" + mol;
					writer.write (pre + ";" + suf + m);
					writer.newLine ();
				}
			}
			if (lexqueFileName != null) writer.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
		ending (fsmFilename);
	}
	
	public void ending (String storeFilename) {
		Enumeration <String> keys = prefixes.keys ();
//		System.out.println ("\n" + prefixes.size () + " states created ");
		idleState.close ("");
//		int count = 0;
		while (keys.hasMoreElements ()) {
			String pref = (String) keys.nextElement ();
			prefixes.get (pref).close (pref);
//			if (count % 1000 == 0) System.out.print ("+");
//			if (count % 100000 == 0) System.out.println (" " + count);
//			count ++;
		}
//		System.out.println ();
		store (storeFilename);
	}
	
	public ArrayList <RecognizedPattern> mark (String src, boolean removeOverlap) {
		FSMState current = prefixes.get ("0");
		if (current == null) current = prefixes.get ("");
		ArrayList <RecognizedPattern> patterns = new ArrayList <RecognizedPattern> ();
//		int tested = 0;
//		int positive = 0;
		for (int i = 0; i < src.length (); i++) {
			char c = src.charAt (i);
			current = current.next (c);
			if (current == null) throw new RuntimeException ("No transition for " + c + ", you should rebuild the FSM including " + c + " into the alphabet");
			if (current.recognized != null) {
				ArrayList <String> suffixes = lexique.get (current.recognized);
				for (String suffixe: suffixes) {
//					tested ++;
					int end = i + 1 + suffixe.length ();
					if (src.length () >= end && src.substring (i + 1, end).equals (suffixe)) {
//						positive ++;
						ArrayList <String> mols = molecules.get (current.recognized + suffixe);
						for (String mol: mols) {
							patterns.add (new RecognizedPattern (mol, i - current.recognized.length () + 1, i + 1 + suffixe.length ()));
						}
					}
				}
			}
		}
//		System.out.println ("Tested: " + tested);
//		System.out.println ("Positive: " + positive);
		if (removeOverlap) {
			ArrayList <RecognizedPattern> patternsFinal = new ArrayList <RecognizedPattern> ();
			for (RecognizedPattern p1: patterns) {
				boolean addToFinal = true;
				for (RecognizedPattern p2: patterns) {
					if (p1 != p2) {
						if (p1.start >= p2.start && p1.start < p2.end || p2.start >= p1.start && p2.start < p1.end) {
							if (p1.length () < p2.length () || p1.length () == p2.length () && p1.start > p2.start) { 
								addToFinal = false;
								break;
							}
						}
					}
				}
				if (addToFinal) {
					patternsFinal.add (p1);
				}
			}
			return patternsFinal;
		}
		else {
			return patterns;
		}
	}
	
	public void store (String filename) {
		if (filename == null) return;
		try {    
            BufferedWriter writer = new BufferedWriter (new FileWriter (filename));
            writer.write (prefixes.size () + "");
			writer.newLine ();
            Enumeration <String> keys = prefixes.keys ();
            keys = prefixes.keys ();
            while (keys.hasMoreElements ()) {
    			String pref = (String) keys.nextElement ();
    			FSMState st = prefixes.get (pref);
    			for (Transition tr: st.transistions) {
    				writer.write (st.id + "�" + tr.read + "�" + tr.nextState.id + "�" + (tr.nextState.recognized == null ? "" : tr.nextState.recognized) + "�");
    			}
    			writer.newLine ();
    		}
            writer.close ();
        } 
        catch (IOException ex) { 
             ex.printStackTrace ();
        } 
	}
	
	public boolean load (String fsmfilename, String lexicFileName) {
		prefixes.clear ();
		molecules.clear ();
		lexique.clear ();
		try {
			File file = new File (fsmfilename);
			if (!file.exists ()) return false;
			
			String l;
			BufferedReader reader = new BufferedReader (new FileReader (lexicFileName));
			while ((l = reader.readLine ()) != null) {
				String ls [] = l.split (";");
				ArrayList <String> mols = new ArrayList <String> ();
				for (int i = 2; i < ls.length; i ++) mols.add (ls [i]);
				molecules.put (ls [0] + ls [1], mols);
				ArrayList <String> sufs = lexique.get (ls [0]);
				if (sufs == null) {
					sufs = new ArrayList <String> ();
					lexique.put (ls [0], sufs);
				}
				sufs.add (ls [1]);
			}
			reader.close ();
			
			reader = new BufferedReader (new FileReader (file));
			int nb = Integer.parseInt (reader.readLine ());
			for (int i = 0; i < nb; i++) {
				prefixes.put (i + "", new FSMState ());
			}
			idleState = prefixes.get ("");
			for (int i = 0; i < nb; i++) {
				l = reader.readLine ();
				String trans [] =  l.split ("�");
				for (int j = 0; j < trans.length; j += 4) {
					FSMState current = prefixes.get (trans [j]);
					FSMState next = prefixes.get (j + 2 >= trans.length ? "" : trans [j + 2]);
					Transition tr = new Transition (trans [j + 1].charAt (0), next);
					if (j + 3 < trans.length && !trans [j + 3].equals ("")) {
						next.recognized = trans [j + 3];
					}
					current.transistions.add (tr);
				}
			}
			reader.close ();
		}
		catch (IOException e) {
			e.printStackTrace ();
		} 
		return true;
	}
	
	/*public static void main (String []args) {
		int max = 18062;
		ArrayList <Integer> lst = new ArrayList <Integer> ();
		for (int i = 0; i < 20; i ++) {
			lst.add ((int) (max * Math.random ())); 
		}
		Collections.sort (lst);
		try {
			BufferedReader reader = new BufferedReader (new FileReader ("workspace/RNA_NCBI_Retrieving_Data/complete_bacterian_genomes.save"));
			int count = 0;
			for (int i = 0; i < lst.size (); i++) {
				String s;
				int v = lst.get (i);
				while ((s = reader.readLine ()) != null) {
					count ++;
					if (count == v) {
						System.out.println (s);
						break;
					}
				}
			}
			
			reader.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}*/
	
//	public static void main (String [] args) {
//		HashFSMBuilder builder = new HashFSMBuilder();
//		builder.normalizeFastaRRNAFile ("workspace/RNA_Test/misc_RNA.fasta");
//		builder.normalizeFastaRRNAFile ("workspace/RNA_Test/mRNA.fasta");
//		builder.normalizeFastaRRNAFile ("workspace/RNA_Test/ncRNA.fasta");
//		builder.normalizeFastaRRNAFile ("workspace/RNA_Test/precursor_RNA.fasta");
//		builder.normalizeFastaRRNAFile ("workspace/RNA_Test/tmRNA.fasta");
//	}
//	
}
