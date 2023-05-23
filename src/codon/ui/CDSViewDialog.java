package codon.ui;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.font.FontRenderContext;
import java.awt.geom.AffineTransform;

import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

import codon.Starter;
import codon.data.ORF;
import codon.recognizer.OrfRecognizer;

public class CDSViewDialog extends JDialog {
	private static final long serialVersionUID = 1L;
	
	JTextField tId = new JTextField ();
	JTextField tEntry = new JTextField ();
	JTextField tName = new JTextField ();
	JTextField tProductId = new JTextField ();
	JTextField tGene = new JTextField ();
	JTextField tIdentity = new JTextField ();
	JComboBox <String> tGo = new JComboBox<String> ();
	JComboBox <String> tKegg = new JComboBox<String> ();
	JTextField tStart = new JTextField ();
	JTextField tStop = new JTextField ();
	JTextField tStartQuery = new JTextField ();
	JTextField tEndQuery = new JTextField ();
	JTextField tStartMatch = new JTextField ();
	JTextField tEndMatch = new JTextField ();
	
//	JTextArea area = new JTextArea ();
	
	AlignmentStringView alignment = new AlignmentStringView ();
	
	boolean viewSequenceAsAmido;
	ORF orf;
	
	class AlignmentStringView extends JPanel {
		private static final long serialVersionUID = 1L;
		String query = null;
		String subject = null;
		Font font = new Font (Font.MONOSPACED, Font.PLAIN, 12);
		
		public AlignmentStringView () {
			JPopupMenu popup = new JPopupMenu ();
			JMenuItem menuCopyQuery = new JMenuItem ("Copy query      Ctrl+Q");
			JMenuItem menuCopySubject = new JMenuItem ("Copy subject   Ctrl+S");
			popup.add (menuCopyQuery);
			popup.add (menuCopySubject);
			setFocusable (true);
			addMouseListener (new MouseAdapter () {
				public void mousePressed (MouseEvent e) {
					grabFocus ();
					if (e.getButton () == MouseEvent.BUTTON3) {
						popup.show (AlignmentStringView.this, e.getX (), e.getY ());
					}
				}
			});
			menuCopyQuery.addActionListener (new ActionListener () {
				public void actionPerformed (ActionEvent e) {
					copyToClipBoard (query);
				}
			});
			menuCopySubject.addActionListener (new ActionListener () {
				public void actionPerformed (ActionEvent e) {
					copyToClipBoard (subject);
				}
			});
			addKeyListener (new KeyAdapter () {
				public void keyPressed (KeyEvent e) {
					if (e.isControlDown () && e.getKeyCode () == KeyEvent.VK_S) {
						copyToClipBoard (subject);
					}
					if (e.isControlDown () && e.getKeyCode () == KeyEvent.VK_Q) {
						copyToClipBoard (query);
					}
				}
			});
		}
				
		public void display (ORF orf) {
			int [] bound = orf.getQueryMatchBounds ();
			int startQuerySeq = bound [0];
			int startMatchSeq = bound [2];
			query = viewSequenceAsAmido ? orf.getAminoSequence () : (orf.isReverse () ? OrfRecognizer.reverseComplementar (orf.getSequence ()) : orf.getSequence ());
			subject = viewSequenceAsAmido ? orf.product.matchingSequence : OrfRecognizer.translateFromAminoToBase (orf.product.matchingSequence);
//			System.err.println (query.substring (startQuerySeq, bound [1]));
//			System.err.println (subject.substring (startMatchSeq, bound [3]));
			
			int startQuery = (startMatchSeq - startQuerySeq) * (viewSequenceAsAmido ? 1 : 3);
			if (startQuery > 0) for (int i = 0; i < startQuery; i++) query = ' ' + query;
			else for (int i = 0; i < -startQuery; i++) subject = ' ' + subject;
			
			FontRenderContext fontRenderer = new FontRenderContext (new AffineTransform (), true, false);
			double qW = font.getStringBounds (query, fontRenderer).getWidth ();
			double sW = font.getStringBounds (subject, fontRenderer).getWidth ();
			setPreferredSize (new Dimension ((int) Math.max (qW, sW) + 10, 80));
			
//			System.err.println ("QUERY AL " + query);
//			System.err.println ("SUBJECT AL " + subject);
		}
		
		public void paintComponent (Graphics g) {
			super.paintComponent (g);
			((Graphics2D) g).setRenderingHint (RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
			g.setFont (font);
			g.drawString (query, 5, 20);
			g.drawString (subject, 5, 40);
		}
	}
	
	protected void copyToClipBoard (String s) {
		 Toolkit.getDefaultToolkit ().getSystemClipboard ().setContents (new StringSelection (s), null);;
	}
	
	public void updateData () {
		
		int acc = (int) (orf.getAccuracy () * 100);
		String title = orf.product.productName + "  (acc: " + acc / 100. + "%)";
		if (orf.product.gene.length () > 0) title = orf.product.gene + " - " + title;
		setTitle (title);
		
		tId.setText (orf.id);
//		tId.setEnabled (false);
		tEntry.setText (orf.product.entryId);
//		tEntry.setEnabled (false);
		tName.setText (orf.product.productName);
		tProductId.setText (orf.product.productId);
//		tProductId.setEnabled (false);
		tGene.setText (orf.product.gene);
//		tGene.setEnabled (false);
		tIdentity.setText (orf.getIdentity () + "");
//		tIdentity.setEnabled (false);
		String s [] = orf.product.goTerm.split ("\\|");
		tGo.removeAllItems ();
		for (int i = 0; i < s.length; i++) {
			tGo.addItem (s [i]);
		}
		if (s.length > 0) tGo.setToolTipText (s [0]);
		s = orf.product.kegg.split ("\\|");
		tKegg.removeAllItems ();
		for (int i = 0; i < s.length; i++) tKegg.addItem (s [i]);
		tStart.setText (orf.start + " (o: " + orf.startOriginal + ")");
//		tStart.setEnabled (false);
		tStop.setText (orf.stop + " (o: " + orf.stopOriginal + ")");
//		tStop.setEnabled (false);
		
		int [] bound = orf.getQueryMatchBounds ();
		int startQuerySeq = bound [0];
		int endQuerySeq = bound [1];
		int startMatchSeq = bound [2];
		int endMatchSeq = bound [3];
		
		tStartQuery.setText (startQuerySeq + " (o: " + orf.product.getStartQuerySeq () + ")");
//		tStartQuery.setEnabled (false);
		tEndQuery.setText (endQuerySeq + " (o: " + orf.product.getEndQuerySeq () + ")");
//		tEndQuery.setEnabled (false);
		tStartMatch.setText (startMatchSeq + " (o: " + orf.product.getStartMatchSeq () + ")");
//		tStartMatch.setEnabled (false);
		tEndMatch.setText (endMatchSeq + " (o: " + orf.product.getEndMatchSeq () + ")");
//		tEndMatch.setEnabled (false);
		
//		String query = viewSequenceAsAmido ? orf.getAmidoSequence () : (orf.isReverse () ? AFDReconhecedorRNA.reverseComplementar (orf.getSequence ()) : orf.getSequence ());
//		String subject = viewSequenceAsAmido ? orf.product.matchingSequence : AFDReconhecedorRNA.translateFromAminoToBase (orf.product.matchingSequence);
//
//		int startQuery = (startMatchSeq - startQuerySeq) * (viewSequenceAsAmido ? 1 : 3);
//		if (startQuery > 0) for (int i = 0; i < startQuery - 1; i++) query = ' ' + query;
//		else for (int i = 0; i < -startQuery; i++) subject = ' ' + subject;
//		
//		area.setText (query + "\n" + subject);
		alignment.display (orf);
	}
	
	public CDSViewDialog (PlotOrfsPanel paneRef, ORF orf, boolean viewSequenceAsAmido) {
		this.orf = orf;
		this.viewSequenceAsAmido = viewSequenceAsAmido;
		setIconImage (Starter.icon);
		orf.loadProduct ();
		
		updateData ();
		
		JPanel pLabelsBase = new JPanel (new GridLayout (6, 1));
		pLabelsBase.setPreferredSize (new Dimension (125, 180));
		pLabelsBase.setMaximumSize (new Dimension (125, 180));
		pLabelsBase.add (new JLabel ("Id"));
		pLabelsBase.add (new JLabel ("Entry"));
		pLabelsBase.add (new JLabel ("Name"));
		pLabelsBase.add (new JLabel ("Product Id"));
		pLabelsBase.add (new JLabel ("Gene"));
		pLabelsBase.add (new JLabel ("Identity"));
		JPanel pFieldsBase = new JPanel (new GridLayout (6, 1));
		pFieldsBase.setPreferredSize (new Dimension (125, 180));
		pFieldsBase.setMaximumSize (new Dimension (125, 180));
		pFieldsBase.add (tId);
		pFieldsBase.add (tEntry);
		pFieldsBase.add (tName);
		pFieldsBase.add (tProductId);
		pFieldsBase.add (tGene);
		pFieldsBase.add (tIdentity);
		
		JPanel pLabelsMatch = new JPanel (new GridLayout (6, 1));
		pLabelsMatch.setPreferredSize (new Dimension (125, 180));
		pLabelsMatch.setMaximumSize (new Dimension (125, 180));
		pLabelsMatch.add (new JLabel ("Start"));
		pLabelsMatch.add (new JLabel ("Stop"));
		pLabelsMatch.add (new JLabel ("startQuery"));
		pLabelsMatch.add (new JLabel ("endQuery"));
		pLabelsMatch.add (new JLabel ("startMatch"));
		pLabelsMatch.add (new JLabel ("endMatch"));
		JPanel pFieldsMatch = new JPanel (new GridLayout (6, 1));
		pFieldsMatch.setPreferredSize (new Dimension (125, 180));
		pFieldsMatch.setMaximumSize (new Dimension (125, 180));
		pFieldsMatch.add (tStart);
		pFieldsMatch.add (tStop);
		pFieldsMatch.add (tStartQuery);
		pFieldsMatch.add (tEndQuery);
		pFieldsMatch.add (tStartMatch);
		pFieldsMatch.add (tEndMatch);
		
		JPanel nPanel = new JPanel (new GridLayout (1, 4));
//		nPanel.setMaximumSize (new Dimension (500, 175));
		nPanel.add (pLabelsBase);
		nPanel.add (pFieldsBase);
		nPanel.add (pLabelsMatch);
		nPanel.add (pFieldsMatch);

		JPanel pTerms = new JPanel (new BorderLayout ());
//		pTerms.setMaximumSize (new Dimension (500, 50));
		JPanel pTermsHead = new JPanel (new GridLayout (2, 1));
		pTermsHead.add (new JLabel ("GO Term"));
		pTermsHead.add (new JLabel ("EC number"));
		JPanel pTermsBody = new JPanel (new GridLayout (2, 1));
		pTermsBody.add (tGo);
		tGo.setPreferredSize (new Dimension (375, 25));
		tGo.setMaximumSize (new Dimension (375, 25));
		tGo.addItemListener (new ItemListener () {
			public void itemStateChanged (ItemEvent e) {
				if (e.getStateChange () == ItemEvent.SELECTED) {
					tGo.setToolTipText (tGo.getSelectedItem ().toString ());
				}
			}
		});
		
		pTermsBody.add (tKegg);
		pTerms.add (pTermsHead, BorderLayout.CENTER);
		pTerms.add (pTermsBody, BorderLayout.EAST);
		
		JPanel mPanel = new JPanel (new BorderLayout ());
		mPanel.add (nPanel);
		mPanel.add (pTerms, BorderLayout.SOUTH);
		
//		area.setPreferredSize (new Dimension (300, 300));
//		area.setColumns (20);
//		area.setRows (5);	
		JScrollPane scroll = new JScrollPane (alignment);
//		area.setEditable (true);
//		area.setFont (new Font (Font.MONOSPACED, Font.PLAIN, 12));
		
		Container root = getContentPane ();
		
		root.add (scroll, BorderLayout.SOUTH);
		if (orf.isBlasted ()) {
			root.add (mPanel, BorderLayout.WEST);
			CDSSelectorPane c = new CDSSelectorPane (paneRef.reco, orf);
			c.addCDSSelectorPaneListener (new CDSSelectorPaneListener () {
				public void productChanged () {
					updateData ();
					paneRef.repaint ();
				}
			});
			JScrollPane sc2 = new JScrollPane (c);
			sc2.setPreferredSize (new Dimension (450, 300));
			root.add (sc2, BorderLayout.CENTER);
			setSize (960, 350);
		}
		else {
			root.add (nPanel);
			setSize (500, 350);
		}
		
		setLocationRelativeTo (paneRef);
		setVisible (true);
	}
}
