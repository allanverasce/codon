package codon.ui.component;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

public class JSelectorMenuItem extends JMenuItem {

	private static final long serialVersionUID = 1L;

	JSpinner spinner = new JSpinner ();
	
	public JSelectorMenuItem (String txt, int threehold) {
		setLayout (new BorderLayout ());
		JLabel label = new JLabel (txt);
		SpinnerModel sm = new SpinnerNumberModel (threehold, 0, 100, 1);
		spinner.setModel (sm);
		add (label);
		add (spinner, BorderLayout.EAST);
		setFocusable (true);
		addMouseListener (new MouseAdapter () {
			public void mousePressed (MouseEvent arg0) {
				grabFocus ();
			}
		});
		setPreferredSize (new Dimension (200, 25));
	}
	
	public int getValue () {
		return (int) spinner.getValue ();
	}
}
