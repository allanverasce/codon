package codon.ui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import codon.data.Codon;
import codon.data.ORF;

public class SearchPanel extends JPanel {
	private static final long serialVersionUID = 1L;
	JTextField tSearch = new JTextField ();
	JLabel lSearch = new JLabel ("Search for:");
	JCheckBox ckEntries = new JCheckBox ("Search in entries");
	JCheckBox ckProducts = new JCheckBox (" Search in other products");
	JCheckBox ckGenes = new JCheckBox (" Search in genes");
	JButton bSearch = new JButton ("Search");
	JList <Integer> lstPositions = new JList <> ();	
	Hashtable <Integer, ORF> searchResults = new Hashtable <Integer, ORF> ();
	PlotOrfsPanel plotOrfs;
	JCheckBox ckExactSentence = new JCheckBox ("Search for the exact sentence");
	JCheckBox ckBase = new JCheckBox ("Base sequence");
	JCheckBox ckAmino = new JCheckBox ("Amino sequence");
	
	JCheckBox ckLowAccuracy = new JCheckBox ("Accuracy <");
	JSpinner spLowAccuracy = new JSpinner ();
	
	JCheckBox ckIntegenicRegion = new JCheckBox ("Interg. Regions > ");
	JSpinner spIntegenicRegion = new JSpinner ();
	
	ActionListener deselectBA =  new ActionListener () {
		public void actionPerformed (ActionEvent arg0) {
			ckBase.setSelected (false);
			ckAmino.setSelected (false);
			ckLowAccuracy.setSelected (false);
			ckIntegenicRegion.setSelected (false);
		}
	};
	
	KeyListener enterPressed = new KeyAdapter () {
		public void keyPressed (KeyEvent e) {
			if (e.getKeyCode () == KeyEvent.VK_ENTER) {
				search ();
			}
		}
	};
	
	public SearchPanel (PlotOrfsPanel plotOrfs) {
		setLayout (new BorderLayout ());
		this.plotOrfs = plotOrfs;
		ckGenes.setSelected (true);
		
		JPanel pOpcoes = new JPanel ();
		pOpcoes.setLayout (new BoxLayout (pOpcoes, BoxLayout.PAGE_AXIS));
		pOpcoes.add (lSearch);
		pOpcoes.add (tSearch);

		JPanel pEntries = new JPanel (new GridLayout (4, 1));
		pEntries.setBorder (BorderFactory.createTitledBorder ("Search in entries"));
		pEntries.add (ckExactSentence);
		pEntries.add (ckGenes);
		pEntries.add (ckProducts);
		pEntries.add (ckEntries);
		ckGenes.addKeyListener (enterPressed);
		ckProducts.addKeyListener (enterPressed);
		ckEntries.addKeyListener (enterPressed);
		ckExactSentence.addKeyListener (enterPressed);
		pOpcoes.add (pEntries);
		
		JPanel pSequence = new JPanel (new GridLayout (2, 1));
		pSequence.setBorder (BorderFactory.createTitledBorder ("Search in sequence"));
		pSequence.add (ckBase);
		pSequence.add (ckAmino);
		ckBase.addKeyListener (enterPressed);
		ckAmino.addKeyListener (enterPressed);
		pOpcoes.add (pSequence);
		
		JPanel pOther = new JPanel (new GridLayout (2, 1));
		pOther.setBorder (BorderFactory.createTitledBorder ("Other criteria"));
		JPanel pLowAccuracy = new JPanel ();
		pLowAccuracy.setLayout (new BoxLayout (pLowAccuracy, BoxLayout.LINE_AXIS));
		pLowAccuracy.add (ckLowAccuracy);
		ckLowAccuracy.addKeyListener (enterPressed);
		pLowAccuracy.add (Box.createHorizontalGlue ());
		pLowAccuracy.add (spLowAccuracy);
		spLowAccuracy.setPreferredSize (new Dimension (70, 22));
		spLowAccuracy.setMaximumSize (new Dimension (70, 22));
		SpinnerModel sm = new SpinnerNumberModel (90, 0, 100, 1);
		spLowAccuracy.setModel (sm);
		spLowAccuracy.addChangeListener (new ChangeListener () {
			public void stateChanged (ChangeEvent arg0) {
				search ();
			}
		});
		JPanel pIntegenicRegion = new JPanel ();
		pIntegenicRegion.setLayout (new BoxLayout (pIntegenicRegion, BoxLayout.LINE_AXIS));
		pIntegenicRegion.add (ckIntegenicRegion);
		ckIntegenicRegion.addKeyListener (enterPressed);
		pIntegenicRegion.add (Box.createHorizontalGlue ());
		pIntegenicRegion.add (spIntegenicRegion);
		spIntegenicRegion.addChangeListener (new ChangeListener () {
			public void stateChanged (ChangeEvent arg0) {
				search ();
			}
		});		
		spIntegenicRegion.setPreferredSize (new Dimension (70, 22));
		sm = new SpinnerNumberModel (500, 0, 100000, 1);
		spIntegenicRegion.setModel (sm);
		pOther.add (pLowAccuracy);
		pOther.add (pIntegenicRegion);
		pOpcoes.add (pOther);
		
		
		ckEntries.addActionListener (deselectBA);
		ckGenes.addActionListener (deselectBA);
		ckProducts.addActionListener (deselectBA);
	
		ckBase.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				ckGenes.setSelected (false);
				ckProducts.setSelected (false);
				ckEntries.setSelected (false);
				ckAmino.setSelected (false);
				ckLowAccuracy.setSelected (false);
				ckIntegenicRegion.setSelected (false);
			}
		});
		
		ckAmino.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				ckGenes.setSelected (false);
				ckProducts.setSelected (false);
				ckEntries.setSelected (false);
				ckBase.setSelected (false);
				ckLowAccuracy.setSelected (false);
				ckIntegenicRegion.setSelected (false);
			}
		});
		
		ckLowAccuracy.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				ckGenes.setSelected (false);
				ckProducts.setSelected (false);
				ckEntries.setSelected (false);
				ckBase.setSelected (false);
				ckBase.setSelected (false);
				ckAmino.setSelected (false);
				ckIntegenicRegion.setSelected (false);
				search ();
			}
		});
		
		ckIntegenicRegion.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				ckGenes.setSelected (false);
				ckProducts.setSelected (false);
				ckEntries.setSelected (false);
				ckBase.setSelected (false);
				ckBase.setSelected (false);
				ckAmino.setSelected (false);
				ckLowAccuracy.setSelected (false);
				search ();
			}
		});
		
		JScrollPane scroll = new JScrollPane (lstPositions);

		add (pOpcoes, BorderLayout.NORTH);
		add (scroll, BorderLayout.CENTER);
		add (bSearch, BorderLayout.SOUTH);
		
		lstPositions.addListSelectionListener (new ListSelectionListener () {
			public void valueChanged (ListSelectionEvent ev) {
				
			}
		});
		lstPositions.addMouseListener (new MouseAdapter () {
			public void mouseClicked (MouseEvent ev) {
				if (ev.getClickCount () == 2) {
					if (lstPositions.getSelectedIndex () != -1) {
						ORF orf = searchResults.get (lstPositions.getSelectedValue ());
						if (orf == null) {
							plotOrfs.setSelectedPosition (lstPositions.getSelectedValue (), true);
						}
						else {
							plotOrfs.selectOrf (searchResults.get (lstPositions.getSelectedValue ()));							
//							plotOrfs.selectOrf (searchResults.get (lstPositions.getSelectedValue ()));
						}						
					}
				}
			}
		});
		tSearch.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				search ();
			}
		});
		bSearch.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				search ();
			}
		});
		
	}
	
	public void search () {
		String v = tSearch.getText ();
		Vector <Integer> positions= new Vector <> ();
		searchResults.clear ();
		if (ckBase.isSelected ()) {
			if (v.length () < 9) {
				JOptionPane.showMessageDialog (PlotOrfsFrame.instance, "Select a sentene higher than 9 bases");
				return;
			}
			StringBuffer seq = plotOrfs.reco.wholeSequence;
			int pos = -1;
			while ((pos = seq.indexOf (v, pos)) != -1) {
				positions.add (pos);
				pos ++;
			}
		}
		else if (ckAmino.isSelected ()) {
			if (v.length () < 4) {
				JOptionPane.showMessageDialog (PlotOrfsFrame.instance, "Select a sentene higher than 4 amino acides");
				return;
			}
			for (int i = 0; i < plotOrfs.reco.amidoSequence.length; i ++) {
				String seq = plotOrfs.reco.amidoSequence [i];
				int pos = -1;
				while ((pos = seq.indexOf (v, pos)) != -1) {
					positions.add (pos * 3 + i % 3);
					pos ++;
				}
			}
		}
		else if (ckLowAccuracy.isSelected ()) {
			for (int i = 0; i < plotOrfs.reco.orfsPorFrame.length; i++) {
				List <ORF> orfs = plotOrfs.reco.orfsPorFrame [i];
				for (ORF orf: orfs) {
					orf.loadProduct ();
					if (orf.isRemoved ()) continue;
					if (orf.getAccuracy () < (int) spLowAccuracy.getValue ()) {
						positions.add (orf.start);
						searchResults.put (orf.start, orf);
					}
				}
			}
		}
		else if (ckIntegenicRegion.isSelected ()) {
			ArrayList <Codon> inter = plotOrfs.reco.intergenicRegions;
			for (int i = 0; i < inter.size (); i += 2) {
				int lg = inter.get (i + 1).index - inter.get (i).index;
				if (lg > (int) spIntegenicRegion.getValue ()) positions.add (inter.get (i).index);
			}
		}
		else {
			if (v.length () == 0) return;
			for (int i = 0; i < plotOrfs.reco.orfsPorFrame.length; i++) {
				List <ORF> orfs = plotOrfs.reco.orfsPorFrame [i];
				for (ORF orf: orfs) {
					orf.loadProduct ();
					if (orf.product.productId == null || orf.product.entryId == null) continue;
					if (ckExactSentence.isSelected ()) {
						if (ckGenes.isSelected () && orf.product.productId.equals (v) || ckProducts.isSelected () && orf.product.productName.equals (v) || ckEntries.isSelected () && orf.product.entryId.equals (v)) {
							positions.add (orf.start);
							searchResults.put (orf.start, orf);
						}
					}
					else {
						if (ckGenes.isSelected () && orf.product.productId.contains (v) || ckProducts.isSelected () && orf.product.productName.contains (v) || ckEntries.isSelected () && orf.product.entryId.contains (v)) {
							positions.add (orf.start);
							searchResults.put (orf.start, orf);
						}
					}
				}
			}
		}
		lstPositions.setListData (positions);
		if (positions.size () != 0) {
			lstPositions.setSelectedIndex (0);
		}
		bSearch.setText ("Search (" + positions.size () + ")");
	}
}

