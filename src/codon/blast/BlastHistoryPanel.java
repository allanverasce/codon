package codon.blast;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.Timer;

public class BlastHistoryPanel extends JPanel {
	private static final long serialVersionUID = 1L;
	
	ArrayList <Integer> blastPerMinute = new ArrayList <> ();
	int last = 0;
	double max = 0;
	
	public BlastHistoryPanel () {
		setPreferredSize (new Dimension ((blastPerMinute.size () + 1) * 12, 60));
		timer.setRepeats (true);
		timer.start ();
	}
	
	Timer timer = new Timer (60000, new ActionListener () {
		public void actionPerformed (ActionEvent arg0) {
			blastPerMinute.add (last);
			if (last > max) max = last;
			last = 0;
			setPreferredSize (new Dimension ((blastPerMinute.size () + 1) * 12, 60));
			updateUI ();
		}
	});
	
	public void newBlast () {
		last ++;
		if (last > max) max = last;
		repaint ();
	}
	
	protected void paintComponent (Graphics g) {
		super.paintComponent (g);
		int w = getWidth ();
		int h = getHeight ();
		int vh = (int) (h * last / max);
		g.drawString (((int) max) + "", getParent ().getWidth ()  / 2 - 4, 12);
		g.fillRect (1, h - vh, 10, vh);
		for (int i = blastPerMinute.size () - 1; i >= 0; i--) {
			int v = blastPerMinute.get (i);
			vh = (int) (h * v / max);
			g.fillRect ((blastPerMinute.size () - i) * 12 + 1, h - vh, 10, vh);
		}
		
		g.setColor (Color.LIGHT_GRAY);
		for (int i = 0; i < max; i += 5) {
			vh = (int) (h * i / max);
			g.drawLine (0, h - vh, w, h - vh);
		}
	}
	
}
