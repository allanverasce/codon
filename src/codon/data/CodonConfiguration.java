package codon.data;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import codon.Starter;
import codon.ui.PlotOrfsFrame;

public class CodonConfiguration extends JDialog {
	private static final long serialVersionUID = 1L;
	
	public static final boolean ENABLE_SWISS = false;
	
	public static String workspace = "workspace/";
	static {
		new File (workspace).mkdirs ();
	}
	
	public static int tolerance = 2;
	public static int overlap_auto_removed = 20;
	public static int accuracy_level_to_remove_overlap_during_blast = 90;
	public static int accuracy_level_to_remove_blasted_orfs_during_blast = 80;
	public static int min_orf_length = 80;
	public static PrintStream out = null;
	
	public static boolean ui = true;
	
	static ArrayList <CodonConfigurationListener> listeners = new ArrayList <> ();
	public static void addCodonConfigurationListener (CodonConfigurationListener l) {
		listeners.add (l);
	}
	public static void removeCodonConfigurationListener (CodonConfigurationListener l) {
		listeners.remove (l);
	}
	public static void fireMinOrfLengthChanged (int minOrfLength) {
		for (CodonConfigurationListener l: listeners) l.minOrfLengthChanged (minOrfLength);
	}
	
	JSpinner spTolerance = new JSpinner ();
	JSpinner spOverlap = new JSpinner ();
	JSpinner spPromote = new JSpinner ();
	JSpinner spDiscard = new JSpinner ();
	JSpinner spOrfMinLenght = new JSpinner ();

	public CodonConfiguration (PlotOrfsFrame frame) {
		setTitle ("Parameters");
		setIconImage (Starter.icon);
		
		Container root = getContentPane ();
	
		JPanel pLabel = new JPanel (new GridLayout (5, 1));
		JPanel pValue = new JPanel (new GridLayout (5, 1));
		JPanel pInfo = new JPanel (new GridLayout (5, 1));
		pLabel.add (new JLabel ("Minimum ORF length"));
		pLabel.add (new JLabel ("Tolerance between accuracies"));
		pLabel.add (new JLabel ("Threshold for overlap auto-remove"));
		pLabel.add (new JLabel ("Accuracy level to promote ORF"));
		pLabel.add (new JLabel ("Accuracy level to discard ORF"));
		
		SpinnerModel sm = new SpinnerNumberModel (tolerance, 0, 100, 1);
		spTolerance.setModel (sm);
		sm = new SpinnerNumberModel (overlap_auto_removed, 0, 100, 1);
		spOverlap.setModel (sm);
		sm = new SpinnerNumberModel (accuracy_level_to_remove_overlap_during_blast, 0, 100, 1);
		spPromote.setModel (sm);
		sm = new SpinnerNumberModel (accuracy_level_to_remove_blasted_orfs_during_blast, 0, 100, 1);
		spDiscard.setModel (sm);
		sm = new SpinnerNumberModel (min_orf_length, 0, 1000, 1);
		spOrfMinLenght.setModel (sm);
		
		pValue.add (spOrfMinLenght);
		pValue.add (spTolerance);
		pValue.add (spOverlap);
		pValue.add (spPromote);
		pValue.add (spDiscard);

		spOrfMinLenght.addChangeListener (new ChangeListener () {
			public void stateChanged(ChangeEvent e) {
				fireMinOrfLengthChanged ((int) spOrfMinLenght.getValue ());
			}
		});
		
		spOrfMinLenght.setToolTipText ("<html>Orf smaller than this threshold will be discarded</html>");
		spTolerance.setToolTipText ("<html>The parameter is used to compare entries.<BR>"
				+ "When the accuracy difference is lower than this threshold, it is not considered as significative.<BR> "
				+ "Then, for instance, an entry with a lower accuracy than another (with a accuracy difference down to this threshold), <BR>but that macth with a gene when the other do not, will be considered better.</html>");
		spOverlap.setToolTipText ("<html>If the overlap percentage is lower than this threshold, the overlap remotion filter will try to rebound<BR> the ORFs but will not remove completly any ORF.</html>");
		spPromote.setToolTipText ("<html>During the blast, if the accuracy of a blasted ORF is higher than this threshold and if the ORF is characterized, <BR>the algorithm will try to remove the overlaps by rebounding or removeing the other unblasted ORF.</html>");
		spDiscard.setToolTipText ("<html>During the blast, if the accuracy of a blasted ORF is higher than the Promote threshold and if the ORF is characterized, <BR>the algorithm will try to remove the ORF that provoke an overlap and that havae an accuracy lower than this threshold.</html>");
		
		ImageIcon intIco = new ImageIcon ("data/img/intA.png");
		JButton iOrfMinLenght = new JButton (intIco);
		JButton iTolerance = new JButton (intIco);
		JButton iOverlap = new JButton (intIco);
		JButton iPromote = new JButton (intIco);
		JButton iDiscard = new JButton (intIco);
		pInfo.add (iOrfMinLenght);
		pInfo.add (iTolerance);
		pInfo.add (iOverlap);
		pInfo.add (iPromote);
		pInfo.add (iDiscard);
		
		iOrfMinLenght.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				JOptionPane.showMessageDialog (CodonConfiguration.this, "<html>Orf smaller than this threshold will be discarded</html>");
			}
		});
		iTolerance.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				JOptionPane.showMessageDialog (CodonConfiguration.this, "<html>The parameter is used to compare entries.<BR>"
						+ "When the accuracy difference is lower than this threshold, it is not considered as significative.<BR> "
						+ "Then, for instance, an entry with a lower accuracy than another (with a accuracy difference down to this threshold), <BR>but that macth with a gene when the other do not, will be considered better.</html>");
			}
		});
		iOverlap.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				JOptionPane.showMessageDialog (CodonConfiguration.this, "<html>If the overlap percentage is lower than this threshold, the overlap remotion filter will try to rebound<BR> the ORFs but will not remove completly any ORF.</html>");
			}
		});
		iPromote.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				JOptionPane.showMessageDialog (CodonConfiguration.this, "<html>During the blast, if the accuracy of a blasted ORF is higher than this threshold and if the ORF is characterized, <BR>the algorithm will try to remove the overlaps by rebounding or removeing the other unblasted ORF.</html>");
			}
		});
		iDiscard.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				JOptionPane.showMessageDialog (CodonConfiguration.this, "<html>During the blast, if the accuracy of a blasted ORF is higher than the Promote threshold and if the ORF is characterized, <BR>the algorithm will try to remove the ORF that provoke an overlap and that havae an accuracy lower than this threshold.</html>");
			}
		});
		
		JButton bCLose = new JButton ("Save and Close");
		bCLose.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				tolerance = (int) spTolerance.getValue ();
				overlap_auto_removed = (int) spOverlap.getValue ();
				accuracy_level_to_remove_overlap_during_blast = (int) spPromote.getValue ();
				accuracy_level_to_remove_blasted_orfs_during_blast = (int) spDiscard.getValue ();
				setVisible (false);
			}
		});
		
		root.add (pLabel, BorderLayout.WEST);
		root.add (pValue, BorderLayout.CENTER);
		root.add (pInfo, BorderLayout.EAST);
		root.add (bCLose, BorderLayout.SOUTH);
		
		setLocationRelativeTo (frame);
		pack ();
		setDefaultCloseOperation (DISPOSE_ON_CLOSE);
		setModal (true);
		setVisible (true);
	}
	
}
