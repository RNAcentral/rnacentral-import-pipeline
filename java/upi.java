public class Upi {

	public static String getUpi( long id ) {
		String str = Long.toHexString( id ).toUpperCase();
		return "URS0000000000".substring(0, 13 - str.length() ) + str;
	}

	public static void main (String args[]) {
		long l = Long.parseLong(args[0]);
		System.out.println(getUpi(l));
	}
}