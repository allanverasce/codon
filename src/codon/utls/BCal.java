package codon.utls;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.GregorianCalendar;

public class BCal extends GregorianCalendar {
	private static final long serialVersionUID = 1L;
	
	public BCal () {}
	
	public BCal (Date dt) {
		super ();
		setTime (dt);
	}
	
	public void addOpenDays (int days) {
        for (int i = 0; i < days; i++) {
            add (DAY_OF_MONTH, 1);
            if (get (DAY_OF_WEEK) == SATURDAY || get (DAY_OF_WEEK) == SUNDAY) i--;
        }
    }
	
	public int getDay () {
		return get (BCal.DAY_OF_MONTH);
	}
	
	public String getMonth () {
		return Meses.MESES [get (MONTH)];
	}
	public int getNumericMonth () {
		return get (MONTH) + 1;
	}
	public String getShortMonth () {
		return Meses.MESES [get (MONTH)].substring (0, 3) + "." ;
	}
	public String getYear () {
		return get (YEAR) + "";
	}
	public String getShortYear () {
		String y = new String (get (YEAR) + "");
		return y.substring (y.length() - 2);
	}
	public String getSqlData () {
		return sqlDate (getTime ());
	}
	public String getSimpleDate () {
		return simpleDate (getTime ());
	}
	static public String sqlDate (Date date) {
		SimpleDateFormat format = new SimpleDateFormat ("yyyy-MM-dd");
		return format.format (date);
	}
	
	static public String simpleDate (Date date) {
		SimpleDateFormat format = new SimpleDateFormat ("dd/MM/yyyy");
		return format.format (date);
	}
	
	static public String shortMes (Date date) {
		SimpleDateFormat format = new SimpleDateFormat ("MM/yy");
		return format.format (date);
	}
	
	static public String completeDate (Date date) {
		BCal cal = new BCal ();
		cal.setTime (date);
		String s = cal.get (DAY_OF_MONTH) + " de " + cal.getMonth () + " de " + cal.getYear ();
		return s;
	}
	
	static public Date getDateFromText (String d) {
		BCal cal = new BCal ();
		cal.set (DAY_OF_MONTH, Integer.parseInt (d.substring (0, d.indexOf ('/'))));
		d = d.substring (d.indexOf ('/') + 1);
		cal.set (MONTH, Integer.parseInt (d.substring (0, d.indexOf ('/'))) - 1);
		d = d.substring (d.indexOf ('/') + 1);
		cal.set (YEAR, Integer.parseInt (d));
		return cal.getTime ();
	}
	
	public String getDayOfWeek () {
		switch (get (BCal.DAY_OF_WEEK)) {
			case BCal.MONDAY: return "Segunda";
			case BCal.TUESDAY: return "Terça";
			case BCal.WEDNESDAY: return "Quarta";
			case BCal.THURSDAY: return "Quinta";
			case BCal.FRIDAY: return "Sexta";
			case BCal.SATURDAY: return "Sábado";
			default: return "Domingo";
		}
	}
	
}
