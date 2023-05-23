package codon.utls;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ContigToList {
	public static void main (String [] args) {
		try {
			BufferedReader reader = new BufferedReader (new FileReader ("workspace/META_GENOMA.contig"));
			BufferedWriter writer = new BufferedWriter (new FileWriter ("workspace/META_GENOMA.fasta"));
			writer.write ("> META_GENOMA");
			writer.newLine ();
			String line = null;
			while ((line = reader.readLine ()) != null) {
				if (line.charAt (0) != '>') writer.write (line);
			}
			reader.close ();
			writer.close ();
			
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
}
