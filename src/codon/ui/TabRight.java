package codon.ui;

import java.awt.Dimension;

import javax.swing.JTabbedPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import codon.blast.BlastingDialog;
import codon.blast.BlastingDialogMonitored;
import codon.data.ORF;


public class TabRight extends JTabbedPane {
	private static final long serialVersionUID = 1L;
	
	MarksPane markPanel = null;
	BlastingDialog singleBlast;
	public OverlapNavigator overlapPane;
	public OrganismPanel organismPane;

	public TabRight (PlotOrfsFrame plotOrfsFrame) {
		markPanel = new MarksPane (plotOrfsFrame.plotOrfsPanel);
		overlapPane = new OverlapNavigator (plotOrfsFrame.plotOrfsPanel);
		organismPane = new OrganismPanel (plotOrfsFrame);
		addTab ("Overlap", overlapPane);
		addTab ("Search", new SearchPanel (plotOrfsFrame.plotOrfsPanel));
		addTab ("Marks", markPanel);
		singleBlast = new BlastingDialog (plotOrfsFrame, 10, false);
		addTab ("Blast", singleBlast);
		addTab ("Organisms", organismPane);
		setPreferredSize (new Dimension (230, 800));
		addChangeListener (new ChangeListener () {
			public void stateChanged (ChangeEvent arg0) {
				if (getSelectedComponent () == organismPane) {
					organismPane.updateFamillies ();
				}
			}
		});
	}
	
	public void addMark (String comment) {
		markPanel.addMark (comment);
	}

	public void addBlast (ORF orf) {
		singleBlast.addOrfToBlast (orf);
	}
	
	public void addBlastTab (BlastingDialog blastingTab) {
		addTab ("Blast all", blastingTab);
		setSelectedComponent (blastingTab);
	}
	
	public void removeBlastTab (BlastingDialog blastingTab) {
		remove (blastingTab);
	}
	
	public void addBlastTab (BlastingDialogMonitored blastingTab) {
		addTab ("Blast all", blastingTab);
		setSelectedComponent (blastingTab);
	}
	
	public void removeBlastTab (BlastingDialogMonitored blastingTab) {
		remove (blastingTab);
	}
}
