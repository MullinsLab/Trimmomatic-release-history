package org.usadellab.trimmomatic.fastq.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class LeadingTrimmer extends AbstractSingleRecordTrimmer
{
	private int qual;

	public LeadingTrimmer(String args)
	{
		qual=Integer.parseInt(args);
	}

	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		String seq=in.getSequence();
		int quals[]=in.getQualityAsInteger(true);
		
		for(int i=0;i<seq.length();i++)
			{
			if(quals[i]>=qual)
				return new FastqRecord(in,i,seq.length()-i);
			}

		return null;
	}

}
