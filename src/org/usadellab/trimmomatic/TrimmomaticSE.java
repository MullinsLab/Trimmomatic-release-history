package org.usadellab.trimmomatic;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.fastq.FastqSerializer;
import org.usadellab.trimmomatic.fastq.trim.Trimmer;
import org.usadellab.trimmomatic.fastq.trim.TrimmerFactory;


public class TrimmomaticSE
{

	/**
	 * Trimmomatic: The FASTQ trimmer
	 * 
	 * CROP:<LENGTH> Crop the read to specified length, by cutting off the right end
	 * LEADING:<QUAL> Trim the read by cutting off the left end while below specified quality
	 * TRAILING:<QUAL> Trim the read by cutting off the right end while below specified quality
	 * SLIDINGWINDOW:<QUAL>:<COUNT> Trim the read once the total quality of COUNT bases drops below QUAL, then trim trailing bases below QUAL 
	 * 
	 * MINLEN:<LENGTH> Drop the read if less than specified length 
	 */

	
	public TrimmomaticSE()
	{
	
	}


	
	public void process(File input, File output, Trimmer trimmers[], int phredOffset, File trimLog) throws IOException
	{
		FastqParser parser=new FastqParser(phredOffset);
		parser.parse(input);

		FastqSerializer serializer=new FastqSerializer();
		serializer.open(output);
		
		PrintStream trimLogStream=null;
		if(trimLog!=null)
			trimLogStream=new PrintStream(trimLog);
		
		FastqRecord recs[]=new FastqRecord[1];
		FastqRecord originalRecs[]=new FastqRecord[1];
		
		while(parser.hasNext())
			{
			originalRecs[0]=recs[0]=parser.next();
			
			for(int i=0;i<trimmers.length;i++)
				recs=trimmers[i].processRecords(recs);
			
			if(recs[0]!=null)
				{
				serializer.writeRecord(recs[0]);
				}
			
			if(trimLogStream!=null)
				{
				for(int i=0;i<originalRecs.length;i++)
					{
					int length=0;
					int startPos=0;
					int endPos=0;
					int trimTail=0;
					
					if(recs[i]!=null)
						{
						length=recs[i].getSequence().length();
						startPos=recs[i].getHeadPos();
						endPos=length+startPos;
						trimTail=originalRecs[i].getSequence().length()-endPos;
						}
					
					trimLogStream.printf("%s %d %d %d %d\n",originalRecs[i].getName(),length,startPos,endPos,trimTail);
					}
				}
			}
		
		serializer.close();
		
		if(trimLogStream!=null)
			trimLogStream.close();
	}
	


	public static void main(String[] args) throws IOException
	{
		int argIndex=0;
		int phredOffset=64;
				
		boolean badOption=false;
		
		File trimLog=null;

		while(args[argIndex].startsWith("-") && argIndex<args.length)
			{
			String arg=args[argIndex++];
			if(arg.equals("-phred33"))
				phredOffset=33;
			else if(arg.equals("-phred64"))
				phredOffset=64;
			else if(arg.equals("-trimlog"))
				{
				if(argIndex<args.length)
					trimLog=new File(args[argIndex++]);
				else
					badOption=true;
				}
			else 
				{
				System.out.println("Unknown option "+arg);
				badOption=true;
				}
			}
	
		if(args.length-argIndex<3 || badOption)
			{
			System.out.println("Usage: TrimmomaticSE [-phred33|-phred64] [-trimlog <trimLogFile>] <inputFile> <outputFile> <trimmer1>...");
			System.exit(1);
			}
		
		File input=new File(args[argIndex++]);
		File output=new File(args[argIndex++]);
		
		TrimmerFactory fac=new TrimmerFactory();
		Trimmer trimmers[]=new Trimmer[args.length-argIndex];
		
		for(int i=0;i<trimmers.length;i++)
			trimmers[i]=fac.makeTrimmer(args[i+argIndex]);
		
		TrimmomaticSE tm=new TrimmomaticSE();
		tm.process(input,output,trimmers,phredOffset,trimLog);
	}

}
