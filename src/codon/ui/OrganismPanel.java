package codon.ui;

import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.DefaultListCellRenderer;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.ListCellRenderer;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import codon.data.ORF;
import codon.data.Organism;
import codon.data.OrganismFamilly;
import codon.recognizer.OrfRecognizer;

public class OrganismPanel extends JPanel {
	private static final long serialVersionUID = 1L;

	ArrayList <OrganismFamilly> famillies;
	OrfRecognizer reco = null;
	DefaultMutableTreeNode treeRoot = new DefaultMutableTreeNode ("Organisms found");
	JList <ORF> orfList = new JList <> ();
	PlotOrfsFrame frame;
	JTree treeOrganisms;
	
	class MyCellRenderer implements ListCellRenderer <ORF> {
		public Component getListCellRendererComponent (JList<? extends ORF> arg0, ORF arg1, int arg2, boolean arg3,	boolean arg4) {
			if (arg1.product != null) return new DefaultListCellRenderer ().getListCellRendererComponent (arg0, arg1.start + " - " + arg1.product.productName, arg2, arg3, arg4);
			else return new DefaultListCellRenderer ().getListCellRendererComponent (arg0, arg1.start + "", arg2, arg3, arg4);
		}
	}
	
	public OrganismPanel (PlotOrfsFrame frame) {
		setLayout (new GridLayout (2, 1));
		this.frame = frame;
		this.reco = frame.reco;

		treeOrganisms = new JTree (treeRoot);
		JScrollPane treeView = new JScrollPane (treeOrganisms);
		
		orfList.setCellRenderer (new MyCellRenderer ());
		orfList.addMouseListener (new MouseAdapter () {
			public void mouseClicked (MouseEvent me) {
				if (me.getClickCount () == 2) {
					ORF orf = orfList.getSelectedValue ();
					if (orf != null) frame.plotOrfsPanel.selectOrf (orf);
				}
			}
		});
		JScrollPane listView = new JScrollPane (orfList);
		
		treeOrganisms.addMouseListener(new MouseAdapter () {
			public void mousePressed (MouseEvent me) {
				if (me.getClickCount () == 2) {
					TreePath tp = treeOrganisms.getPathForLocation (me.getX (), me.getY ());
					if (tp != null) {
						DefaultMutableTreeNode node = (DefaultMutableTreeNode) tp.getLastPathComponent ();
						Object ob = hash.get(node.toString ());
						if (ob == null) orfList.setListData (new Vector <ORF> ());
						else {
							if (ob instanceof OrganismFamilly) {
								orfList.setListData (((OrganismFamilly) ob).orfs);
							}
							else {
								orfList.setListData (((Organism) ob).orfs);
							}
						}
						orfList.updateUI ();
					}
				}
			}
		});
		
		add (treeView);
		add (listView);
	}
	
	Hashtable <String, Object> hash = new Hashtable<String, Object> (); 
	
	public void updateFamillies () {
		famillies = reco.getOrganisms ();
		treeRoot.removeAllChildren ();
		hash.clear ();
		for (OrganismFamilly fam: famillies) {
			DefaultMutableTreeNode famNode = new DefaultMutableTreeNode (fam.famillyName + " (" + fam.occurrences + " products - " + fam.size () + " organisms)");
			treeRoot.add (famNode);
			hash.put (famNode.toString (), fam);
			for (Organism org: fam) {
				DefaultMutableTreeNode orgNode = new DefaultMutableTreeNode (org.organismName + " (" + org.occurrences + " products)");
				famNode.add (orgNode);
				hash.put (orgNode.toString (), org);
			}
		}
		treeOrganisms.updateUI ();
	}
}
