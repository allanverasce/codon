package codon.ui;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.tree.DefaultMutableTreeNode;

import codon.Starter;
import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.data.Organism;
import codon.data.OrganismFamilly;
import codon.recognizer.OrfRecognizer;

public class ReportDialog extends JDialog {
	private static final long serialVersionUID = 1L;
	
	JLabel lGenes = new JLabel ("Genes: ");
	JLabel lRNAs = new JLabel ("RNA: ");
	JLabel ltRNAs = new JLabel ("tRNA: ");
	JLabel loRNAs = new JLabel ("Other RNA: ");
	JLabel lProducts = new JLabel ("Other products: ");
	JLabel lUncharacterized = new JLabel ("Uncharacterized: ");
	JLabel lTotal = new JLabel ("Total: ");
	JLabel lFree = new JLabel ("Intergenic regions (%): ");
	JLabel lOverlap = new JLabel ("Overlap regions (%): ");
	JLabel lAccuracy = new JLabel ("Accuracy (%): ");
	
	JTextField lvGenes = new JTextField ();
	JTextField lvRNAs = new JTextField ();
	JTextField lvtRNAs = new JTextField ();
	JTextField lvoRNAs = new JTextField ();
	JTextField lvProducts = new JTextField ();
	JTextField lvUncharacterized = new JTextField ();
	JTextField lvTotal = new JTextField ();
	JTextField lvFree = new JTextField ();
	JTextField lvOverlap = new JTextField ();
	JTextField lvAccuracy = new JTextField ();
	
	OrfRecognizer reco;
	
	public ReportDialog (PlotOrfsFrame frame) {
		setIconImage (Starter.icon);
		reco = frame.reco;

		Container root = getContentPane ();
		
		JPanel panelReport = new JPanel (new GridLayout (10, 2));
		
		JMenuBar menuBar = new JMenuBar ();
		setJMenuBar (menuBar);
		JMenu menuFile = new JMenu ("File");
		JMenuItem menuExportBlastResuts = new JMenuItem ("Export blast informations");
		menuBar.add (menuFile);
		menuFile.add (menuExportBlastResuts);
		menuExportBlastResuts.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				JFileChooser chooser = new JFileChooser (CodonConfiguration.workspace);
				FileNameExtensionFilter filter = new FileNameExtensionFilter ("CSV", "csv");
				chooser.setFileSelectionMode (JFileChooser.FILES_ONLY);
				chooser.setFileFilter (filter);
			    if (chooser.showOpenDialog (ReportDialog.this) == JFileChooser.APPROVE_OPTION) {
			    	saveBlastResult (chooser.getSelectedFile ());
			    }
			}
		});
		
		reco.computeOrfFreeAndOrfSuperpositionSpans (0);
		
		double freeSpanLength = 0;
		double overlapSpanLength = 0;
		for (int i = 0; i < reco.intergenicRegions.size (); i += 2) {
			freeSpanLength += reco.intergenicRegions.get (i + 1).index - reco.intergenicRegions.get (i).index - 2;
		}
		for (int i = 0; i < reco.overlapSpans.size (); i += 2) {
			overlapSpanLength += reco.overlapSpans.get (i + 1).index - reco.overlapSpans.get (i).index;
		}
		
		reco.computeOrfFreeAndOrfSuperpositionSpans (frame.tabRight.overlapPane.getTolerance ());
		
		int genes = 0;
		int total = 0;
		int products = 0;
		int uncharacterized = 0;
		double accuracy = 0;
		int rrnas = 0;
		int trnas = 0;
		int ornas = 0;
		for (int i = 0; i < 6; i ++) {
			List <ORF> orfs = reco.orfsPorFrame [i];
			for (ORF orf: orfs) {
				if (!orf.isRemoved ()) {
					if (orf.isGene ()) genes ++;
					else if (orf.product.entryId.equals ("RNA")) rrnas ++;
					else if (orf.product.entryId.equals ("tRNA")) trnas ++;
					else if (orf.product.entryId.contains ("RNA")) ornas ++;
					else if (!orf.isUncharacterized ()) products ++;
					else uncharacterized ++;
					total ++;
					accuracy += orf.getAccuracy ();
				}
			}
		}
		
				
//		products -= reco.rnaOrfs.size ();
		accuracy /= total;
		freeSpanLength /= reco.wholeSequence.length ();
		overlapSpanLength /= reco.wholeSequence.length ();
		lvGenes.setText ("" + genes);
		lvGenes.setEditable (false);
		lvRNAs.setText ("" + rrnas);
		lvRNAs.setEditable (false);
		lvtRNAs.setText ("" + trnas);
		lvtRNAs.setEditable (false);
		lvoRNAs.setText ("" + ornas);
		lvoRNAs.setEditable (false);
		lvProducts.setText ("" + products);
		lvProducts.setEditable (false);
		lvUncharacterized.setText ("" + uncharacterized);
		lvUncharacterized.setEditable (false);
		lvTotal.setText ("" + total);
		lvTotal.setEditable (false);
		lvFree.setText ("" + String.format ("%.4f", freeSpanLength));
		lvFree.setEditable (false);
		lvOverlap.setText ("" + String.format ("%.4f", overlapSpanLength));
		lvOverlap.setEditable (false);
		lvAccuracy.setText ("" + String.format ("%.2f", accuracy));
		lvAccuracy.setEditable (false);
		setTitle ("ORFs: " + total);
		
		panelReport.add (lGenes);
		panelReport.add (lvGenes);
		panelReport.add (lRNAs);
		panelReport.add (lvRNAs);
		panelReport.add (ltRNAs);
		panelReport.add (lvtRNAs);
		panelReport.add (loRNAs);
		panelReport.add (lvoRNAs);
		panelReport.add (lProducts);
		panelReport.add (lvProducts);
		panelReport.add (lUncharacterized);
		panelReport.add (lvUncharacterized);
		panelReport.add (lTotal);
		panelReport.add (lvTotal);
		
		panelReport.add (lFree);
		panelReport.add (lvFree);
		panelReport.add (lOverlap);
		panelReport.add (lvOverlap);
		panelReport.add (lAccuracy);
		panelReport.add (lvAccuracy);
		
		ArrayList <OrganismFamilly> famillies = reco.getOrganisms ();
		DefaultMutableTreeNode treeRoot = new DefaultMutableTreeNode ("Organisms found");
		for (OrganismFamilly fam: famillies) {
			DefaultMutableTreeNode famNode = new DefaultMutableTreeNode (fam.famillyName + " (" + fam.occurrences + " products - " + fam.size () + " organisms)");
			treeRoot.add (famNode);
			for (Organism org: fam) {
				DefaultMutableTreeNode orgNode = new DefaultMutableTreeNode (org.organismName + " (" + org.occurrences + " products)");
				famNode.add (orgNode);
			}
		}
		JTree treeOrganisms = new JTree (treeRoot);
		JScrollPane treeView = new JScrollPane (treeOrganisms);
		
		root.add (panelReport);
		root.add (treeView, BorderLayout.SOUTH);
		
		pack ();
		setLocationRelativeTo (frame);
		setVisible (true);
	}

	protected void saveBlastResult (File f) {
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter (f));
			writer.write ("ORF id;UNIPROT Entry ID; Initial length;Final length;Blasted;Saved blasts;Is gene;Is uncharacterized protein;Accuracy");
			writer.newLine ();
			for (ArrayList <ORF> orfs: reco.orfsPorFrame) for (ORF orf: orfs) {
				if (orf.length > CodonConfiguration.min_orf_length 
						|| (orf.stopOriginal - orf.startOriginal) > CodonConfiguration.min_orf_length
						|| orf.product != null && orf.product.entryId != null && orf.product.entryId.contains ("RNA")) {
					writer.write (	orf.id + ";" +
									(orf.product != null && orf.product.entryId != null? orf.product.entryId : "") + ";" +
									orf.isRemoved () + ";" +
									(orf.stopOriginal - orf.startOriginal) + ";" + 
									(orf.stop - orf.start) + ";" + 
									orf.isBlasted () + ";" +
									orf.killed.size () + ";" +
									orf.isGene () + ";" +
									orf.isUncharacterized () + ";" +
									String.format ("%2f", orf.getAccuracy ()));
					writer.newLine ();
				}
			}
			writer.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
}
