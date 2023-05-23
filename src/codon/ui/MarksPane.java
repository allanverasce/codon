package codon.ui;

import java.awt.BorderLayout;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;

import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import codon.Starter;
import codon.data.CodonConfiguration;

public class MarksPane extends JPanel {
	private static final long serialVersionUID = 1L;
	JList <String> list = new JList <String> ();
	Vector <String> data = new Vector <> ();
	
	PlotOrfsPanel plotOrfsPanel;
	public MarksPane (PlotOrfsPanel plotOrfsPanel) {
		this.plotOrfsPanel = plotOrfsPanel;
		try {
			File f = new File (CodonConfiguration.workspace + plotOrfsPanel.reco.dataDestDirS + "/marked.txt");
			if (f.exists ()) {
				BufferedReader reader = new BufferedReader (new FileReader (f));
				String s = null;
				while ((s = reader.readLine ()) != null) {
					data.add (s);
				}
				reader.close ();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
		list.setListData (data);
		JScrollPane scroll = new JScrollPane (list);
		setLayout (new BorderLayout ());
		add (scroll);
		list.addKeyListener (new KeyAdapter () {
			public void keyPressed (KeyEvent ev) {
				if (ev.getKeyCode () == KeyEvent.VK_DELETE && list.getSelectedIndex () != -1) {
					data.remove (list.getSelectedIndex ());
					list.setListData (data);
					try {
						BufferedWriter writer = new BufferedWriter(new FileWriter (CodonConfiguration.workspace + plotOrfsPanel.reco.dataDestDirS + "/marked.txt"));
						for (int i = 0; i < list.getModel ().getSize (); i++) {
							writer.write (list.getModel ().getElementAt (i));
							writer.newLine ();
						}
						writer.close ();
					} catch (IOException e) {
						e.printStackTrace ();
					}
				}
			}
		});
		list.addMouseListener (new MouseAdapter () {
			public void mouseClicked (MouseEvent arg0) {
				if (list.getSelectedIndex () != -1) {
					String s = list.getSelectedValue ();
					s = s.substring (0, s.indexOf (':'));
					int pos = Integer.parseInt (s);
					plotOrfsPanel.setSelectedPosition (pos, true);
				}
			}
		});
	}
	
	public void addMark (String comment) {
		data.add (comment);
		list.updateUI ();
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter (CodonConfiguration.workspace + plotOrfsPanel.reco.dataDestDirS + "/marked.txt", true));
			writer.write (comment);
			writer.newLine ();
			writer.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
}
