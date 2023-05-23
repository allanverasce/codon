package codon.ui;

import java.awt.BorderLayout;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JDialog;
import javax.swing.JLabel;

import org.apache.commons.io.FileUtils;

import codon.Starter;
import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.recognizer.OrfRecognizer;

class ImportORFSplash extends JDialog {
	private static final long serialVersionUID = 1L;
	
	JLabel lSearchFor = new JLabel ("Search for 10000/10000");
	JLabel lFound = new JLabel ("Found 0");
	int searchFor = 0;
	int found = 0;
	
	public ImportORFSplash (PlotOrfsFrame frame, OrfRecognizer reco, OrfRecognizer reference) {
		getContentPane ().add (lSearchFor, BorderLayout.CENTER);
		getContentPane ().add (lFound, BorderLayout.SOUTH);
		setUndecorated (true);
		pack ();
		setLocationRelativeTo (frame);
		setModal (true);
		
		for (int i = 0; i < reference.orfsPorFrame.length; i++) {
			ArrayList <ORF> orfs = reference.orfsPorFrame [i];
			for (ORF orf: orfs) {
				if (!orf.isRemoved ()) {
					searchFor ++;
				}
			}
		}
		lSearchFor.setText ("Search for 1/" + searchFor);
		
		Thread th = new Thread () {
			public void run () {
				int indice = 0;
				for (int f = 0; f < reference.orfsPorFrame.length; f++) {
					ArrayList <ORF> orfs = reference.orfsPorFrame [f];
					for (ORF orf: orfs) {
						if (!orf.isRemoved ()) {
							indice ++;
							lSearchFor.setText ("Search for " + indice + "/" + searchFor);
							String seq = orf.getSequence ();
							int last = 0;
							while ((last = reco.wholeSequence.indexOf (seq, last)) != -1) {
								found ++;
								lFound.setText ("Found " + found);
								int frame = (last - 2) % 3 + (f > 2 ? 3 : 0);
								ArrayList <ORF> existingOrfs = reco.orfsPorFrame [frame];
								orf.loadProduct ();
								for (int j = 0; j < existingOrfs.size (); j++) {
									ORF eorf = existingOrfs.get (j);
									if (frame < 3 && (eorf.stop + 3 == (last + seq.length ()))) {
										eorf.product = orf.product;
										eorf.setLoaded (true);
										eorf.start = eorf.stop + 3 - seq.length ();
										try {
											FileUtils.copyFile (new File (CodonConfiguration.workspace + reference.dataDestDirS + "/orfs/" + orf.id + ".unip"), new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + eorf.id + ".unip"));
											if (orf.hasSwissUpdate ()) {
												FileUtils.copyFile (new File (CodonConfiguration.workspace + reference.dataDestDirS + "/orfs/" + orf.id + ".swiss"), new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + eorf.id + ".swiss"));
											}
										} catch (IOException e) {
											e.printStackTrace ();
										}
										break;
									}
									if (frame > 2 && eorf.start == last) {
										eorf.product = orf.product;
										eorf.setLoaded (true);
										eorf.stop = eorf.start + seq.length ();
										try {
											FileUtils.copyFile (new File (CodonConfiguration.workspace + reference.dataDestDirS + "/orfs/" + orf.id + ".unip"), new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + eorf.id + ".unip"));
											if (orf.hasSwissUpdate ()) {
												FileUtils.copyFile (new File (CodonConfiguration.workspace + reference.dataDestDirS + "/orfs/" + orf.id + ".swiss"), new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + eorf.id + ".swiss"));
											}
										} catch (IOException e) {
											e.printStackTrace ();
										}
										break;
									}
									
								}
								last += 1;
							}
						}
					}
				}
				frame.updateFilter ();
				setVisible (false);
			}
		};	
		th.start ();
		
		setVisible (true);
	}
	
	
}