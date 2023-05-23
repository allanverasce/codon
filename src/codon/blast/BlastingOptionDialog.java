package codon.blast;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.TitledBorder;

import codon.Starter;
import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.recognizer.OrfRecognizer;
import codon.ui.PlotOrfsFrame;

public class BlastingOptionDialog extends JDialog {
	private static final long serialVersionUID = 1L;
	
	JLabel lStartPos = new JLabel ("First sequence position");
	JTextField spStartPos = new JTextField ();
	JLabel lEndPos = new JLabel ("Last sequence position");
	JTextField spEndPos = new JTextField ();
	JLabel lThreadPerFrame = new JLabel ("Simultanous requests");
//	JTextField spThreadPerFrame = new JTextField ();
	JSpinner spThreadPerFrame = new JSpinner ();
	
	JButton bStart = new JButton ("Start blasting");
	ArrayList <ORF> validOrfs = new ArrayList <> ();
	
	JCheckBox ckAccelerate = new JCheckBox ("Accelerate with SwissProt");
	
	OrfRecognizer reco;
	boolean fullBlast;
	
	class PanelBlasted extends JPanel {
		private static final long serialVersionUID = 1L;
		public PanelBlasted () {
			setPreferredSize (new Dimension (validOrfs.size (), 30));
		}
		public void paint (Graphics g) {
			int count = 0;
			for (int i = 0; i < validOrfs.size (); i++) {
				if (validOrfs.get (i).isBlasted ()) g.setColor (Color.WHITE);
				else {
					g.setColor (Color.RED);
					count ++;
				}
				g.drawLine (i, 0, i, 30);
			}
			g.setColor (Color.BLACK);
			for (int i = 0; i < validOrfs.size (); i += 500) {
				g.drawLine (i, 0, i, 30);
				g.drawString (i + "", i + 5, 20);
			}
			g.drawString (count + "", 100, 20);
		}
	}
	
	public BlastingOptionDialog (PlotOrfsFrame plotOrfsFrame, OrfRecognizer reco, boolean fullBlast) {
		setTitle ("Blasting options");
		setIconImage (Starter.icon);
		this.reco = reco;
		this.fullBlast = fullBlast;
		for (int frame = 0; frame < 6; frame++) {
			List <ORF> orfs = reco.orfsPorFrame [frame];
			for (ORF orf: orfs) {
				if (!orf.isRemoved ()) validOrfs.add (orf);
			}
		}
		
		JPanel pOptions = new JPanel (new GridLayout (7, 1));
		TitledBorder title = BorderFactory.createTitledBorder ("Options");
		pOptions.setBorder (title);
		pOptions.add (lStartPos);
		pOptions.add (spStartPos);
		spStartPos.setText ("0");
		pOptions.add (lEndPos);
		pOptions.add (spEndPos);
		spEndPos.setText (reco.wholeSequence.length () + "");
		pOptions.add (lThreadPerFrame);
		pOptions.add (spThreadPerFrame);
		spThreadPerFrame.setValue (1);
		if (CodonConfiguration.ENABLE_SWISS) pOptions.add (ckAccelerate);
		ckAccelerate.setSelected (false);
		ckAccelerate.setToolTipText ("<HTML>When selected, the software will try to blast with SwissProt database first.<BR>"
				+ "If there is no hit with high accuracy, it will blast with UniprotKB database.<BR>"
				+ "The blast with SwissProt is faster than with UniprotKB, however more genes are deposited into Uniprot.<BR>"
				+ "The relevance of this option may depend on: the organism studied, internet bandwidth, Uniprot and SwissProt workloads.<BR>"
				+ "The efficacy is unpredictable, whereas the option will be AUTOMATICALLY desactivated if during the last 60 tries,<BR>"
				+ "the algorithm observes no effective time gain.");
		Container root = getContentPane ();
		JScrollPane scrollBlasted = new JScrollPane (new PanelBlasted ());
		scrollBlasted.setPreferredSize (new Dimension (300, 50));
		root.add (scrollBlasted, BorderLayout.NORTH);
		root.add (pOptions, BorderLayout.CENTER);
		root.add (bStart, BorderLayout.SOUTH);
		
		SpinnerModel sm = new SpinnerNumberModel (1, 1, 20, 1);
		spThreadPerFrame.setModel (sm);
		
		bStart.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				int first = 0;
				int last = 0;
				int threadNumber = 0;
				try {
					first = Integer.parseInt (spStartPos.getText ());
				} catch (Exception e) {
					JOptionPane.showMessageDialog (BlastingOptionDialog.this, "First sequence position must be a number");
					return;
				}
				try {
					last = Integer.parseInt (spEndPos.getText ());
				} catch (Exception e) {
					JOptionPane.showMessageDialog (BlastingOptionDialog.this, "Last sequence position must be a number");
					return;
				}
				try {
					threadNumber = (Integer) spThreadPerFrame.getValue ();
				} catch (Exception e) {
					JOptionPane.showMessageDialog (BlastingOptionDialog.this, "Simultanous requests/frame");
					return;
				}
				Vector <ORF> subValid = new Vector <> ();
				for (ORF orf: validOrfs) {
					if (orf.start >= first && orf.start < last) subValid.add (orf);
				}
				new BlastingDialogMonitored (plotOrfsFrame, subValid, threadNumber, ckAccelerate.isSelected (), fullBlast);
				setVisible (false);
			}
		});
		pack ();
		setLocationRelativeTo (plotOrfsFrame);
		setVisible (true);
	}
	
}
