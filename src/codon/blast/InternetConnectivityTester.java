package codon.blast;

import java.io.IOException;

import javax.swing.JOptionPane;

public class InternetConnectivityTester {
	public static boolean isConnected () {
		try {
			String os = System.getProperty ("os.name");
			Process process;
			if (os.contains ("Windows")) {
				process = java.lang.Runtime.getRuntime ().exec ("ping 8.8.8.8");
			}
			else {
				process = java.lang.Runtime.getRuntime ().exec ("ping -c 4 8.8.8.8");
			}
	        int x = process.waitFor ();
			return x == 0;
		} catch (InterruptedException e) {
			e.printStackTrace ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
		return false;
	}
	
	public static void main (String [] args) {
		JOptionPane.showMessageDialog (null, InternetConnectivityTester.isConnected () + "");
	}
	
}
