package demo;

import java.io.File;

import codon.data.CodonConfiguration;
import codon.recognizer.OrfRecognizer;

public class ExpeRNA {
	public static void main (String[] args) {
		CodonConfiguration.ui = false;
		new File ("workspace/RNA_EXPE_RES/").mkdirs ();
		File sourcesFiles [] = new File ("workspace/RNA_EXPE").listFiles ();
		for (File srcF: sourcesFiles) {
			if (srcF.isFile ()) {
				OrfRecognizer reco = new OrfRecognizer ();
				reco.loadFromEmbl (srcF);
				reco.annotateRNAOrf (true);
				reco.exportAsEmbl (new File ("workspace/RNA_EXPE_RES/" + srcF.getName ()));
			}
		}
	}
}
