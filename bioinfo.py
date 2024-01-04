import streamlit as st
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from Bio.Seq import Seq 
from Bio import SeqIO
from collections import Counter
import neatbio.sequtils as utils
import numpy as np
from PIL import Image 
import io


def delta(x,y):
    return 0 if x == y else 1


def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]
def plotMatrix(M,t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)


def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2) #experiment with character choice
# Convert to Fxn
def dotplotx(seq1,seq2):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    xt=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # on y-axis list all sequences of seq 1
    yt=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()


def gc_content(seq):
    result = float(str(seq).count('G') + str(seq).count('C'))/len(seq) * 100
    return result

def at_content(seq):
    result = float(str(seq).count('A') + str(seq).count('T'))/len(seq) * 100
    return result



def main():
    """A Simple Streamlit App """
    st.title("BioInformatics App By Pradeep Ganap⚛⏩")

    activity = ['Intro', 'DNA Sequence', 'DotPlot']
    choice = st.sidebar.selectbox("Select Activity", activity)
    if choice == 'Intro':
        st.subheader("Intro to BioInformatics")
        image = Image.open('dna-logo.jpg')
        
        st.image(image, use_column_width=True)
        st.write("""
        # Simple BioInformatics Web App On The Principle Of Central Dogma
        
        This app helps in 'DNA Sequence Analysis' from the uploded FASTA file, which gives 
        detailed info about Nucleotide Frequency with the PLot Freq, DNA Composition with 
        GC and AT Content, as well as Protein Synthesis including Transcription, Translation,
        Complement Sequence, AA Frequency, Plot for AA Frequency and Full Amino Acid Name.
        
        ## BioInformatics : 
        It is the science concerned with the development and applcation of computer hardware 
        and software to the acquisition,storage,analysis, and visualization of biological information.
                         
        ## Central Dogma : 
        The central dogma of molecular biology is a theory stating that genetic information 
        flows only in one direction, from DNA, to RNA, to Protein, or RNA directly to Protein.
        
        ## Three Processes : 
        From existing DNA to make new DNA (DNA Replcation) 
        From DNA to make new RNA (Transcription)
        From RNA to make new Proteins or Amino Acid Chain (Translation).
        
        ***
        """)
    elif choice == "DNA Sequence":
        st.subheader("DNA Sequence Analysis")
        seq_file = st.file_uploader("Upload FASTA File",type=["fasta","fa"])
        if seq_file is not None:
            byte_str = seq_file.read()
            text_obj = byte_str.decode('UTF-8')
            dna_record = SeqIO.read(io.StringIO(text_obj),"fasta")
            st.write(dna_record)
            # st.write(dna_record)
            dna_seq = dna_record.seq

            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record.description)
            elif details == "Sequence":
                st.write(dna_record.seq)
            # Nucleotide Frequencies
            st.subheader("Nucleotide Frequency")
            dna_freq = Counter(dna_seq)
            st.write(dna_freq)
            adenine_color = st.color_picker("Adenine Color")
            thymine_color = st.color_picker("Thymine Color")
            guanine_color = st.color_picker("Guanine Color")
            cytosine_color = st.color_picker("Cytosine Color")

            if st.button("Plot Freq"):
                barlist = plt.bar(dna_freq.keys(),dna_freq.values())
                barlist[2].set_color(adenine_color)
                barlist[3].set_color(thymine_color)
                barlist[1].set_color(guanine_color)
                barlist[0].set_color(cytosine_color)
                st.pyplot()
                st.set_option('deprecation.showPyplotGlobalUse', False)

            st.subheader("DNA Composition")
            gc_score = utils.gc_content(str(dna_seq))
            at_score = utils.at_content(str(dna_seq))
            st.json({"GC Content":gc_score,"AT Content":at_score})

            # Nucleotide Count
            nt_count = st.text_input("Enter Nucleotide Here","Type Nucleotide Alphabet")
            st.write("Number of {} Nucleotide is ::{}".format((nt_count),str(dna_seq).count(nt_count)))

            # Protein Synthesis
            st.subheader("Protein Synthesis")
            p1 = dna_seq.translate()
            aa_freq = Counter(str(p1))

            if st.checkbox("Transcription"):
                st.write(dna_seq.transcribe())

            elif st.checkbox("Translation"):
                st.write(dna_seq.translate())

            elif st.checkbox("Complement"):
                st.write(dna_seq.complement())

            elif st.checkbox("AA Frequency"):
                st.write(aa_freq)

            elif st.checkbox("Plot AA Frequency"):
                aa_color = st.color_picker("Pick An Amino Acid Color")
                # barlist = plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                # barlist[2].set_color(aa_color)
                plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                st.pyplot()
                fig, ax = plt.subplots()
                ax.scatter([1, 2, 3], [1, 2, 3])
                st.pyplot(fig)
                st.set_option('deprecation.showPyplotGlobalUse', False)

            elif st.checkbox("Full Amino Acid Name"):
                aa_name = str(p1).replace("*","")
                aa3 = utils.convert_1to3(aa_name)
                st.write(aa_name)
                st.write("=====================")
                st.write(aa3)

                st.write("=====================")
                st.write(utils.get_acid_name(aa3))

            # Top Most Common Amino
    elif choice == "DotPlot":
        st.subheader("Generate Dot Plot For Two Sequences")
        seq_file1 = st.file_uploader("Upload 1st FASTA File",type=["fasta","fa"])
        seq_file2 = st.file_uploader("Upload 2nd FASTA File",type=["fasta","fa"])

        if seq_file1 and seq_file2 is not None:
            dna_record1 = SeqIO.read(seq_file1,"fasta")
            dna_record2 = SeqIO.read(seq_file2,"fasta")
            # st.write(dna_record)
            dna_seq1 = dna_record1.seq
            dna_seq2 = dna_record2.seq

            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record1.description)
                st.write("=====================")
                st.write(dna_record2.description)
            elif details == "Sequence":
                st.write(dna_record1.seq)
                st.write("=====================")
                st.write(dna_record2.seq)

            cus_limit = st.number_input("Select Max number of Nucleotide",10,200,50)
            if st.button("Dot Plot"):
                st.write("Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                dotplotx(dna_seq1[0:cus_limit],dna_seq2[0:cus_limit])

                st.pyplot()
        
if __name__ == '__main__':
    main()
