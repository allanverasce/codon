package codon.ui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Vector;

import javax.swing.DefaultListCellRenderer;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import codon.data.Codon;
import codon.recognizer.OrfRecognizer;
import codon.recognizer.OrfRecognizerListener;

public class OverlapNavigator extends JPanel implements OrfRecognizerListener {
	private static final long serialVersionUID = 1L;
	
	DefaultListCellRenderer r;
	
	Vector <Integer> startsArea = new Vector <> ();
	Vector <Boolean> solved = new Vector <> ();
	JList <Integer> lstPositions = new JList <> ();
	JLabel lQuantity = new JLabel ("0/0");
	
	JLabel lTolerance = new JLabel ("Overlap tolerance (bases)");
	JSpinner spiTolerance = new JSpinner ();
	
	public class MyRenderer extends DefaultListCellRenderer {
		private static final long serialVersionUID = 1L;
		public Component getListCellRendererComponent (JList <?> arg0, Object arg1, int index, boolean arg3, boolean arg4) {
			if (arg0 == null) System.err.println ("OverlapNavigator arg0 null");
			if (arg1 == null) System.err.println ("OverlapNavigator arg1 null");
			Component comp = super.getListCellRendererComponent (arg0, arg1, index, arg3, arg4);
			if (index < solved.size ()) return comp;
			if (solved.get (index) == true) comp.setBackground (Color.LIGHT_GRAY);
			return comp;
		}
	}
	
	PlotOrfsPanel plotOrfs;
	public OverlapNavigator (PlotOrfsPanel plotOrfs) {
		setLayout (new BorderLayout ());
		this.plotOrfs = plotOrfs;
		JButton bPrevious = new JButton (new ImageIcon ("./data/img/prevA.png"));
		JButton bNext = new JButton (new ImageIcon ("./data/img/nextA.png"));
		JButton bRefresh = new JButton (new ImageIcon ("./data/img/reloadA.png"));
		
		SpinnerModel sm = new SpinnerNumberModel (10, 0, 100, 1);
		spiTolerance.setModel (sm);
		
		JPanel navPane = new JPanel (new FlowLayout ());
		navPane.add (bPrevious);
		bPrevious.setMargin (new Insets (0, 0, 0, 0));
		bPrevious.setPreferredSize (new Dimension (40, 40));
		navPane.add (bRefresh);
		bRefresh.setMargin (new Insets (0, 0, 0, 0));
		bRefresh.setPreferredSize (new Dimension (45, 40));
		navPane.add (bNext);
		bNext.setMargin (new Insets (0, 0, 0, 0));
		bNext.setPreferredSize (new Dimension (40, 40));
		navPane.add (lQuantity);
		JPanel tolerancePane = new JPanel (new FlowLayout ());
		tolerancePane.add (lTolerance);
		tolerancePane.add (spiTolerance);
		JPanel headerPane = new JPanel (new GridLayout (2, 1));
		headerPane.add (navPane);
		headerPane.add (tolerancePane);
		
		plotOrfs.reco.computeOrfFreeAndOrfSuperpositionSpans ((int) spiTolerance.getValue ());
		refresh ();
		plotOrfs.reco.addOrfRecognizerListener (this);
		
		JScrollPane scroll = new JScrollPane (lstPositions);
		
		lstPositions.addMouseListener (new MouseAdapter () {
			public void mouseClicked (MouseEvent arg0) {
				if (lstPositions.getSelectedIndex () != -1) {
					plotOrfs.setSelectedPosition (lstPositions.getSelectedValue (), true);
				}
			}
		});
		bPrevious.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				if (lstPositions.getSelectedIndex () > 0) {
					int i = lstPositions.getSelectedIndex () - 1;
					while (i != 0 &&  solved.get (i) == true) {
						i --;
					}
					lstPositions.setSelectedValue (startsArea.get (i), true);
				}
				else {
					lstPositions.setSelectedValue (startsArea.get (0), true);
				}
				plotOrfs.setSelectedPosition (lstPositions.getSelectedValue (), true);
			}
		});
		bNext.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				if (lstPositions.getSelectedIndex () != -1 && lstPositions.getSelectedIndex () < startsArea.size () - 1) {
					int i = lstPositions.getSelectedIndex () + 1;
					while (i < startsArea.size () - 1 &&  solved.get (i) == true) {
						i ++;
					}
					lstPositions.setSelectedValue (startsArea.get (i), true);
				}
				else {
					lstPositions.setSelectedValue (startsArea.get (0), true);
				}
				plotOrfs.setSelectedPosition (lstPositions.getSelectedValue (), true);
			}
		});
		bRefresh.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				refresh ();
			}
		});
		spiTolerance.addChangeListener (new ChangeListener () {
			public void stateChanged (ChangeEvent arg0) {
				refresh ();
			}
		});
		
		add (headerPane, BorderLayout.NORTH);
		add (scroll, BorderLayout.CENTER);
		lstPositions.setCellRenderer (new MyRenderer ());
	}
	
	public int getTolerance () {
		return (int) spiTolerance.getValue ();
	}
	
	public void refresh () {
		solved.clear ();
		startsArea.clear ();
		plotOrfs.reco.computeOrfFreeAndOrfSuperpositionSpans ((int) spiTolerance.getValue ());
		for (Codon codon: plotOrfs.reco.overlapSpans) {
			if (codon.type == OrfRecognizer.START) {
				solved.add (false);
				startsArea.add (codon.index);
			}
		}
		lstPositions.setListData (startsArea);
		lQuantity.setText ("0/" + startsArea.size ());
	}

	public void superpositionsAndFreeSpansUpdate () {
		int i = 0;
		int j = 0;
		ArrayList <Codon> lst = plotOrfs.reco.overlapSpans;
		if (lst.size () == startsArea.size () / 2) return;
		solved.clear ();
		int solvedCount = 0;
		while (i < startsArea.size () && j < lst.size ()) {
			if (startsArea.get (i) < lst.get (j).index) {
				solved.add (true);
				solvedCount ++;
				i ++;
			}
			else if (startsArea.get (i) == lst.get (j).index) {
				solved.add (false);
				i ++;
				j += 2;
			}
			else {
				startsArea.add (i, lst.get (j).index);
				solved.add (false);
				j += 2;
			}
		}
		if (solved.size () < startsArea.size ()) for (int k = solved.size (); k < startsArea.size (); k ++) solved.add (true);
		else if (startsArea.size () < lst.size () / 2) for (int k = startsArea.size (); k < lst.size () / 2 - 1; k += 2) {
			startsArea.add (lst.get (k).index);
			solved.add (false);
		}
		lQuantity.setText (solvedCount + "/" + startsArea.size ());
		lstPositions.updateUI ();
	}
	
}
