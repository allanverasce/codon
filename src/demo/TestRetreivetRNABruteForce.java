package demo;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import codon.data.CodonConfiguration;
import codon.recognizer.OrfRecognizer;

public class TestRetreivetRNABruteForce {

	public static void main (String[] args) {
		OrfRecognizer reco = new OrfRecognizer ();
		reco.load (CodonConfiguration.workspace + "PALESB58.fasta");
		String seq = reco.wholeSequence.toString ();
		try {
			long init = System.currentTimeMillis ();
			BufferedReader reader = new BufferedReader (new FileReader ("data/bd/rrna.fasta"));
			String l;
			int i = 0;
			while ((l = reader.readLine ()) != null) {
				l = reader.readLine ();
				int ind = seq.indexOf (l);
				while (ind != -1) ind = seq.indexOf (l, ind + 1);
				if (i % 1000 == 0) System.out.print ("-");
				if (i % 10000 == 0) System.out.println ();
				i ++;
			}
			reader.close ();
			System.out.println ((System.currentTimeMillis () - init) / 1000);
		} catch (FileNotFoundException e) {
			e.printStackTrace ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
}
