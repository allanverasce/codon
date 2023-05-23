package demo;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class CountOrganisms {
	public static void main (String [] args) throws IOException { //14302 //5305
		BufferedReader reader = new BufferedReader (new FileReader ("workspace/RNA_NCBI_Retrieving_Data/all_bacterian_genomes.save.txt"));
		String line;
		Hashtable <String, String> orgs = new Hashtable <String, String> ();
		while ((line = reader.readLine ()) != null) {
			if (line.indexOf (";") != -1) line = line.substring (0, line.indexOf (";"));
			orgs.put (line, line);
		}
		System.out.println (orgs.values ().size ());
		reader.close ();
	}
}
