package codon;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Enumeration;
import java.util.Hashtable;

public class UpgradeRNA {
	
	Hashtable <String, String> rrnas = new Hashtable <String, String> ();
	Hashtable <String, String> trnas = new Hashtable <String, String> ();
	
	public void loadTRNA () {
		try {
			BufferedReader reader = new BufferedReader (new FileReader ("data/bd/trna-bacterial.fasta"));
			String s;
			String molecula = "";
			StringBuffer rna = null;
			int redundant = 0;
			while ((s = reader.readLine ()) != null) {
				if (s.charAt (0) == '>') {
					if (rna != null) {
						String seq = rna.toString ();
						trnas.put (seq, molecula);
						if (trnas.get (seq) != null) redundant ++;
						trnas.put (seq, molecula);
					}
					molecula = s.substring (1, s.indexOf (" "));
					rna = new StringBuffer ();
				}
				else {
					rna.append (s);
				}
			}
			String seq = rna.toString ();
			trnas.put (seq, molecula);
			reader.close ();
			//if (redundant != 0) System.out.println (redundant + " TRNA redundant occurencies ignored");
		}
		catch (Exception e) {
		}
	}
	
	public void loadRRNA () {
		try {
			BufferedReader reader = new BufferedReader (new FileReader ("data/bd/rrna.fasta"));
			String s;
			String molecula = "";
			StringBuffer rna = null;
			int redundant = 0;
			while ((s = reader.readLine ()) != null) {
				if (s.charAt (0) == '>') {
					if (rna != null) {
						String seq = rna.toString ();
						if (rrnas.get (seq) != null) redundant ++;
						rrnas.put (seq, molecula);
					}
					if (s.contains ("/molecule=")) {
						molecula = s.substring (s.indexOf ("=") + 1, s.indexOf ("/s") - 1);
					}
					else {
						molecula = s.substring (1, s.indexOf (" "));
					}
					rna = new StringBuffer ();
				}
				else {
					rna.append (s);
				}
			}
			String seq = rna.toString ();
			rrnas.put (seq, molecula);
			reader.close ();
			//if (redundant != 0) System.out.println (redundant + " RRNA redundant occurencies ignored");
		}
		catch (Exception e) {
		}
		
	}
	
	public UpgradeRNA () {
		loadRRNA ();
		loadTRNA ();
		int initTRNAs = trnas.size ();
		int initRRNAs = rrnas.size ();
		int trnasReceived = 0;
		int trnasAdded = 0;
		int rrnasReceived = 0;
		int rrnasAdded = 0;
		BufferedReader reader = new BufferedReader (new InputStreamReader (System.in));
		String seq;
		try {
			while (!(seq = reader.readLine ()).equals ("##")) {
				String name = reader.readLine ();
				String exists = trnas.get (seq);
				trnasReceived ++;
				if (exists == null) {
					trnasAdded ++;
					//System.out.println ("\t\tadding tRNA:" + name + "  " + seq);
					trnas.put (seq, name);
				}
			}
			while (!(seq = reader.readLine ()).equals ("##")) {
				String name = reader.readLine ();
				String exists = rrnas.get (seq);
				rrnasReceived ++;
				if (exists == null) {
					//System.out.println ("\t\tadding rRNA:" + name + "  " + seq);
					rrnasAdded ++;
					rrnas.put (seq, name);
				}
			}
		}
		catch (IOException io) {
		}
		int endTRNAs = trnas.size ();
		int endRRNAs = rrnas.size ();
		if (initTRNAs != endTRNAs) {
			try {
				BufferedWriter writer = new BufferedWriter (new FileWriter ("data/bd/trna-bacterial.fasta"));
				Enumeration <String> seqs = trnas.keys ();
				while (seqs.hasMoreElements ()) {
					seq = (String) seqs.nextElement ();
					if (seq.length () == 0) continue; 
					String name = trnas.get (seq);
					writer.write (">" + name + " ");
					writer.newLine ();
					writer.write (seq);
					writer.newLine ();
				}
				writer.close ();
			} catch (IOException e) {
				e.printStackTrace ();
			}
			new File ("data/bd/trna-bacterial.fsm").delete ();
		}
		System.out.println ("\ttRNA: begin=" + initTRNAs + " received=" + trnasReceived + " added=" + trnasAdded + " final=" + endTRNAs);
		if (initRRNAs != endRRNAs) {
			try {
				BufferedWriter writer = new BufferedWriter (new FileWriter ("data/bd/rrna.fasta"));
				Enumeration <String> seqs = rrnas.keys ();
				while (seqs.hasMoreElements ()) {
					seq = (String) seqs.nextElement ();
					if (seq.length () == 0) continue; 
					String name = rrnas.get (seq);
					writer.write (">" + name + " ");
					writer.newLine ();
					writer.write (seq);
					writer.newLine ();
				}
				writer.close ();
			} catch (IOException e) {
				e.printStackTrace ();
			}
		}
		System.out.println ("\trRNA: begin=" + initRRNAs + " received=" + rrnasReceived + " added=" + rrnasAdded + " final=" + endRRNAs);
	}
	
	public static void main (String [] args) {
		new UpgradeRNA ();
	}

}
