package agkanogiannis;

import agkanogiannis.d2stools.Read;
import agkanogiannis.d2stools.ReadD2;

public class Coder {
	
	//private static Alphabet DNA = Alphabet.DNA;
	
	public static long encodeToLong(Read read, int start, int end, boolean reverse){
		if(!reverse){	
			long num = 0L;
			for(int i=start; i<end; i++){
				num = num * 4L;
				switch (read.charAt(i)){
					case 'A':
					case 'a':
						num = num + 0L;
					break;
	
					case 'T':
					case 't':
						num = num + 1L;
					break;
	
					case 'C':
					case 'c':
						num = num + 2L;
					break;
	
					case 'G':
					case 'g':
						num = num + 3L;
					break;
	
					default:
						num = -1L;
					break;
				}
				if (num == -1L){
					break;
				}
			}
			return num;
		}
		else{
			long num = 0L;
			for(int i=end-1; i>=start; i--){
				num = num * 4L;
				switch (read.charAt(i)){
					case 'A':
					case 'a':
						num = num + 1L;
					break;
	
					case 'T':
					case 't':
						num = num + 0L;
					break;
	
					case 'C':
					case 'c':
						num = num + 3L;
					break;
	
					case 'G':
					case 'g':
						num = num + 2L;
					break;
	
					default:
						num = -1L;
					break;
				}
				if (num == -1L){
					break;
				}
			}
			return num;
		}
	}
	
	public static long encodeToLong(ReadD2 read, int start, int end, boolean reverse, boolean countNucl){
		if(!reverse){	
			long num = 0L;
			short a=0;
			short t=0;
			short c=0;
			short g=0;
			for(int i=start; i<end; i++){
				num = num * 4L;
				switch (read.charAt(i)){
					case 'A':
					case 'a':
						num = num + 0L;
						if(countNucl)
							a++;
					break;
	
					case 'T':
					case 't':
						num = num + 1L;
						if(countNucl)
							t++;
					break;
	
					case 'C':
					case 'c':
						num = num + 2L;
						if(countNucl)
							c++;
					break;
	
					case 'G':
					case 'g':
						num = num + 3L;
						if(countNucl)
							g++;
					break;
	
					default:
						num = -1L;
					break;
				}
				if (num == -1L){
					break;
				}
			}
			if (countNucl && num != -1L){
				read.as += a;
				read.ts += t;
				read.cs += c;
				read.gs += g;
			}
			return num;
		}
		else{
			long num = 0L;
			short a=0;
			short t=0;
			short c=0;
			short g=0;
			for(int i=end-1; i>=start; i--){
				num = num * 4L;
				switch (read.charAt(i)){
					case 'A':
					case 'a':
						num = num + 1L;
						if(countNucl)
							t++;
					break;
	
					case 'T':
					case 't':
						num = num + 0L;
						if(countNucl)
							a++;
					break;
	
					case 'C':
					case 'c':
						num = num + 3L;
						if(countNucl)
							g++;
					break;
	
					case 'G':
					case 'g':
						num = num + 2L;
						if(countNucl)
							c++;
					break;
	
					default:
						num = -1L;
					break;
				}
				if (num == -1L){
					break;
				}
			}
			if (countNucl && num != -1L){
				read.as += a;
				read.ts += t;
				read.cs += c;
				read.gs += g;
			}
			return num;
		}
	}
	
	public static String decodeFromLong(long encoded, int length){
		StringBuilder sbf = new StringBuilder();
		StringBuilder sbr = new StringBuilder();
		for(int i=0; i<length; i++){
			int mod = (int)(encoded%4L);
			switch(mod){
				case 0:
					sbf.append('A');
					sbr.append('T');
					encoded -= 0L;
					break;
				case 1:	
					sbf.append('T');
					sbr.append('A');
					encoded -= 1L;
					break;
				case 2:	
					sbf.append('C');
					sbr.append('G');
					encoded -= 2L;
					break;
				case 3:	
					sbf.append('G');
					sbr.append('C');
					encoded -= 3L;
					break;
				default:	
					break;
			}
			encoded /= 4L;
		}
		String f =  sbf.reverse().toString();
		String r = sbr.toString();
		return f+"/"+r;
	}
	
	public static long encodeToLong(String kmerS, int k){
		long num = 0L;
		for(int i=0; i<k; i++){
			num = num * 4L;
			switch (kmerS.charAt(i)){
				case 'A':
				case 'a':
					num = num + 0L;
				break;

				case 'T':
				case 't':
					num = num + 1L;
				break;

				case 'C':
				case 'c':
					num = num + 2L;
				break;

				case 'G':
				case 'g':
					num = num + 3L;
				break;

				default:
					num = -1L;
				break;
			}
			if (num == -1L){
				break;
			}
		}
		return num;
	}
	
	/*
	public static byte[] encode(Read read, int start, int end){
		int size = ((end-start-1)>>2) + 1;
		byte[] encoded = new byte[size];
		int mod = (end-start)%4;
		int d1,d2,d3,d4;
		byte b;
		for(int i=0; i<size; i++){
			if(i < size-1){
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				d2 = DNA.toIndex(read.charAt((i<<2) + 1 + start));
				d3 = DNA.toIndex(read.charAt((i<<2) + 2 + start));
				d4 = DNA.toIndex(read.charAt((i<<2) + 3 + start));
				b =  (byte) (d1 + (d2<<2) + (d3<<4) + (d4<<6));
				encoded[i] = b;
				continue;
			}
			switch (mod) {
			case 0:
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				d2 = DNA.toIndex(read.charAt((i<<2) + 1 + start));
				d3 = DNA.toIndex(read.charAt((i<<2) + 2 + start));
				d4 = DNA.toIndex(read.charAt((i<<2) + 3 + start));
				b =  (byte) (d1 + (d2<<2) + (d3<<4) + (d4<<6));
				encoded[i] = b;
				break;
			case 1:
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				b = (byte)d1;
				encoded[i] = b;
				break;
			case 2:
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				d2 = DNA.toIndex(read.charAt((i<<2) + 1 + start));
				b =  (byte) (d1 + (d2<<2));
				encoded[i] = b;
				break;
			case 3:
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				d2 = DNA.toIndex(read.charAt((i<<2) + 1 + start));
				d3 = DNA.toIndex(read.charAt((i<<2) + 2 + start));
				b =  (byte) (d1 + (d2<<2) + (d3<<4));
				encoded[i] = b;
				break;
			default:
				break;
			}
		}
		return encoded;
	}
	*/
	
	/*
	public static char[] decode(byte[] encoded, int k) {
		char[] decoded = new char[k];
		int size = encoded.length;
		int mod = k % 4;
		int d1, d2, d3, d4;
		byte b;
		for (int i = 0; i < size; i++) {
			b = encoded[i];
			if (i < size - 1) {
				d1 = (int) (b & (3));
				d2 = (int) ((b & (3 << 2)) >> 2);
				d3 = (int) ((b & (3 << 4)) >> 4);
				d4 = (int) ((b & (3 << 6)) >> 6);
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				decoded[(i << 2) + 1] = DNA.toChar(d2);
				decoded[(i << 2) + 2] = DNA.toChar(d3);
				decoded[(i << 2) + 3] = DNA.toChar(d4);
				continue;
			}
			switch (mod) {
			case 0:
				d1 = (int) (b & (3));
				d2 = (int) ((b & (3 << 2)) >> 2);
				d3 = (int) ((b & (3 << 4)) >> 4);
				d4 = (int) ((b & (3 << 6)) >> 6);
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				decoded[(i << 2) + 1] = DNA.toChar(d2);
				decoded[(i << 2) + 2] = DNA.toChar(d3);
				decoded[(i << 2) + 3] = DNA.toChar(d4);
				break;
			case 1:
				d1 = (int) (b & (3));
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				break;
			case 2:
				d1 = (int) (b & (3));
				d2 = (int) ((b & (3 << 2)) >> 2);
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				decoded[(i << 2) + 1] = DNA.toChar(d2);
				break;
			case 3:
				d1 = (int) (b & (3));
				d2 = (int) ((b & (3 << 2)) >> 2);
				d3 = (int) ((b & (3 << 4)) >> 4);
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				decoded[(i << 2) + 1] = DNA.toChar(d2);
				decoded[(i << 2) + 2] = DNA.toChar(d3);
				break;
			default:
				break;
			}
		}
		return decoded;
	}
	*/

}