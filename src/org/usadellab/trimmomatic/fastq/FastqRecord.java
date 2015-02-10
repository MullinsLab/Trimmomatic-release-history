package org.usadellab.trimmomatic.fastq;

public class FastqRecord
{
	private String name;
	private String sequence;
	private String comment;
	private String quality;
	
	private int phredOffset;
	private int headPos;
	
	public FastqRecord(String name, String sequence, String comment, String quality, int phredOffset)
	{
		this.name=name;
		this.sequence=sequence;
		this.comment=comment;
		this.quality=quality;
	
		this.phredOffset=phredOffset;
		headPos=0;
		
		if(sequence.length()!=quality.length())
			throw new RuntimeException("Sequence and quality length don't match: '"+sequence+"' vs '"+quality+"'");
	}

	public FastqRecord(FastqRecord base, int headPos, int length)
	{
		this.sequence=base.sequence.substring(headPos,headPos+length);
		this.quality=base.quality.substring(headPos,headPos+length);
		
		this.name=base.name;
		this.comment=base.comment;
		this.phredOffset=base.phredOffset;	
		
		this.headPos=base.headPos+headPos;
	}

	public FastqRecord(FastqRecord base, String sequence, String quality, int phredOffset)
	{
		this.sequence=sequence;
		this.quality=quality;
		
		this.name=base.name;
		this.comment=base.comment;
		this.headPos=base.headPos;

		this.phredOffset=phredOffset;
	}

	public String getName()
	{
		return name;
	}

	public String getSequence()
	{
		return sequence;
	}
	
	public String getComment()
	{
		return comment;
	}

	public String getQuality()
	{
		return quality;
	}
	
	public int getPhredOffset()
	{
		return phredOffset;
	}
	
	public int getHeadPos()
	{
		return headPos;
	}
	

	public int[] getQualityAsInteger(boolean zeroNs)
	{
		int arr[]=new int[quality.length()];
		
		for(int i=0;i<quality.length();i++)
			{
			if(zeroNs && sequence.charAt(i)=='N')
				arr[i]=0;
			else
				arr[i]=quality.charAt(i)-phredOffset;
			}
			
		return arr;
	}
	
	
}
