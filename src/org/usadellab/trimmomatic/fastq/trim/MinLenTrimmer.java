package org.usadellab.trimmomatic.fastq.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class MinLenTrimmer extends AbstractSingleRecordTrimmer
{
	private int minLen;

	public MinLenTrimmer(String args)
	{
		minLen=Integer.parseInt(args);
	}

	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		if(in.getSequence().length()>=minLen)
			return in;
		
		return null;
	}

}
