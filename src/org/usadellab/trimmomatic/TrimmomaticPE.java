package org.usadellab.trimmomatic;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.fastq.FastqSerializer;
import org.usadellab.trimmomatic.fastq.trim.Trimmer;
import org.usadellab.trimmomatic.fastq.trim.TrimmerFactory;


public class TrimmomaticPE
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

	
	public TrimmomaticPE()
	{
	
	}


	
	public void process(File input1, File input2, File output1P, File output1U, File output2P, File output2U, Trimmer trimmers[], int phredOffset, File trimLog) throws IOException
	{
		FastqParser parser1=new FastqParser(phredOffset);
		parser1.parse(input1);

		FastqParser parser2=new FastqParser(phredOffset);		
		parser2.parse(input2);
		
		FastqSerializer serializer1P=new FastqSerializer();
		serializer1P.open(output1P);
		
		FastqSerializer serializer1U=new FastqSerializer();
		serializer1U.open(output1U);
		
		FastqSerializer serializer2P=new FastqSerializer();
		serializer2P.open(output2P);
		
		FastqSerializer serializer2U=new FastqSerializer();
		serializer2U.open(output2U);
		
		FastqRecord originalRecs[]=new FastqRecord[2];
		FastqRecord recs[]=new FastqRecord[2];
		
		PrintStream trimLogStream=null;
		if(trimLog!=null)
			trimLogStream=new PrintStream(trimLog);
		
		while(parser1.hasNext() && parser2.hasNext())
			{
			originalRecs[0]=recs[0]=parser1.next();
			originalRecs[1]=recs[1]=parser2.next();
			
			for(int i=0;i<trimmers.length;i++)
				{
				recs=trimmers[i].processRecords(recs);
				}
			
			if(recs[0]!=null && recs[1]!=null)
				{
				serializer1P.writeRecord(recs[0]);
				serializer2P.writeRecord(recs[1]);
				}
			else if(recs[0]!=null)
				{
				serializer1U.writeRecord(recs[0]);
				}
			else if(recs[1]!=null)
				{
				serializer2U.writeRecord(recs[1]);
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
		
		serializer1P.close();
		serializer1U.close();
		serializer2P.close();
		serializer2U.close();
		
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
	
		if(args.length-argIndex<7 || badOption)
			{
			System.out.println("Usage: TrimmomaticPE [-phred33|-phred64] [-trimlog <trimLogFile>] <inputFile1> <inputFile2> <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U> <trimmer1>...");
			System.exit(1);
			}
		
		File input1=new File(args[argIndex++]);
		File input2=new File(args[argIndex++]);
		
		File output1P=new File(args[argIndex++]);
		File output1U=new File(args[argIndex++]);
		
		File output2P=new File(args[argIndex++]);
		File output2U=new File(args[argIndex++]);
		
		TrimmerFactory fac=new TrimmerFactory();
		Trimmer trimmers[]=new Trimmer[args.length-argIndex];
		
		for(int i=0;i<trimmers.length;i++)
			trimmers[i]=fac.makeTrimmer(args[i+argIndex]);
		
		TrimmomaticPE tm=new TrimmomaticPE();
		tm.process(input1,input2,output1P,output1U,output2P,output2U,trimmers,phredOffset,trimLog);
	}

}
