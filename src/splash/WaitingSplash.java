package splash;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;

import codon.Starter;
import codon.ui.PlotOrfsFrame;

public class WaitingSplash extends JDialog {
	private static final long serialVersionUID = 1L;

	JProgressBar progress;
	JLabel lState;
	
	public WaitingSplash (Component relativeTo, String waitingPlease, String inicialState, int progressBarUpLimit, WaitingSplashThread th, boolean canCancel) {
		JLabel lWaitingPlease = new JLabel (waitingPlease);
		lState = new JLabel (inicialState);
		progress = new JProgressBar (0, progressBarUpLimit);
		JButton bCancel = new JButton ("Cancel");
		
		JPanel p = new JPanel (new GridLayout (3, 1));
		p.add (lWaitingPlease);
		p.add (lState);
		p.add (progress);
		
		getContentPane ().add (p);
		if (canCancel) {
			getContentPane ().add (bCancel, BorderLayout.SOUTH);
			bCancel.addActionListener (new ActionListener () {
				public void actionPerformed (ActionEvent e) {
					th.mustStop = true;
					th.interrupt ();
					setVisible (false);	
				}
			});
		}
		
		getContentPane ().add (new JLabel (new ImageIcon (Starter.iconBordrered)), BorderLayout.NORTH);
		getContentPane ().add (new JPanel (), BorderLayout.WEST);
		getContentPane ().add (new JPanel (), BorderLayout.EAST);
		setUndecorated (true);
		setLocationRelativeTo (relativeTo == null ? PlotOrfsFrame.instance : relativeTo);
		pack ();
		setDefaultCloseOperation (JDialog.DISPOSE_ON_CLOSE);
		setModal (true);
	}
	
	public void changeStatus (String state, int progressValue) {
		SwingUtilities.invokeLater (new Runnable () {
			public void run () {
				lState.setText (state);
				progress.setValue (progressValue);
			}
		});
	}
	
//	static WaitingSplash splash;
//	public static void main (String [] args) {
//		System.err.println ("Activethread at start " + Thread.activeCount ());
//		WaitingSplashThread th = new WaitingSplashThread () {
//			public void run () {
//				for (int i = 0; i < 10000; i++) {
//					try {
//						Thread.sleep (10);
//						if (splash != null) splash.changeStatus (i + "/ 1 zilião", i);
//					} catch (InterruptedException e) {
//						break;
//					}
//				}
//			}
//			
//		};
//		splash = new WaitingSplash (null, "Teste de paciencia", "0/ 1 zilião", 10000, th, true);
//		th.start ();
//		System.err.println ("Activethread now " + Thread.activeCount ());
//		splash.setVisible (true);
//		System.err.println ("Activethread now " + Thread.activeCount ());
//	} 
}
