import os
import streamlit as st
import streamlit.components.v1 as components
from streamlit_hgb import hgb, reference_hash, load_samples, hgb_run
import pandas as pd
import gffutils
import glob
from streamlit_drawable_canvas import st_canvas
import time

DB = "hg38.genes.db"
_RELEASE = True
# app: `$ streamlit run main.py`

# Retrieve gene annottions from gff file
#@st.cache(allow_output_mutation=True)
def load_db(db_file):
    return gffutils.FeatureDB(db_file, keep_order=True)

if _RELEASE:
    st.header("Hybrid Genome Browser")

    # Create a second instance of our component whose `name` arg will vary
    # based on a text_input widget.
    #
    # We use the special "key" argument to assign a fixed identity to this
    # component instance. By default, when a component's arguments change,
    # it is considered a new instance and will be re-mounted on the frontend
    # and lose its current state. In this case, we want to vary the component's
    # "name" argument without having it get recreated.
    # name_input = st.text_input("Enter a file name", value="../../bt142/ont2_ngmlr.bam")
    try:
        yaml = load_samples("config.yaml")

        ref = st.sidebar.selectbox("Which references to use?", list(yaml.keys()), 1)

        name_input = st.sidebar.multiselect("Which files to load?",
          yaml[ref]["alignments"],
          list(yaml[ref]["default"])
        )
        refs = reference_hash(yaml[ref]["alignments"][0]) 
        default_range = yaml[ref]["range"][0]
        db_file = yaml[ref]["db"][0]
    except:
        files = glob.glob("*.bam")
        if len(files) > 0:
            name_input = st.sidebar.multiselect("Which files to load?",
                files,
                [files[0]])
            refs = reference_hash(name_input[0]) 
        else:
            name_input = st.sidebar.text_input("Which file to explore?")
            refs = reference_hash(name_input)
            name_input = [name_input]
        db_file = DB
        if len(refs) > 0:
            #default_range = "{}:10001-20001".format(next(iter(refs)))
            #default_range = "{}:8794744-8850896".format(next(iter(refs)))
            default_range = "{}:8874744-8950896".format(next(iter(refs)))
        else:
            default_range = ""
    
    region = st.sidebar.text_input("Where to explore?", default_range)
    bed = st.sidebar.text_input("Which bed file to use for ?", "")
    split=False
    coverage=20
    y=16
    callet=True
    no_ins=False
    db = load_db(db_file)

    if len(refs) > 0:
        try:
            # Fetch from gene
            gene = db[region] #load_db(db_file, region)
            #print(gene, gene.seqid)
            chr_def = gene.seqid
            car, cdr = gene.start, gene.end
            default_range = "{}:{}-{}".format(chr_def, car, cdr)
        except: 
            try:
                chr_def, region_chr = region.split(":")
                car, cdr = region_chr.split("-")
            except:
                region = default_range
                chr_def, region_chr = region.split(":")
                car, cdr = region_chr.split("-")
        chr = list(refs.keys())

        ref_id = st.sidebar.selectbox(
        'Which chromosome to display?',
        chr, chr.index(chr_def))
        range = st.sidebar.slider(
         'Select a range of values',
         0, refs[ref_id], (int(car), int(cdr)))

        if st.sidebar.checkbox("Detail"):
            num = st.sidebar.number_input("Enter a start coordinate", 0, refs[ref_id], range[0])
            num2 = st.sidebar.number_input("Enter a stop coordinate", 0, refs[ref_id], range[1])
            coverage = st.sidebar.number_input('The expected coverage', 1, 500, coverage)
            split = st.sidebar.checkbox('Split-alignment only view')
            callet = st.sidebar.checkbox('Show callets only intra-chromosomal split alignment', True)
            y = st.sidebar.number_input("Set a read height", 8, 128, y)

        if range[1] - range[0] <= 1000*1000*12:
            if bed != "":
               flags = " -J {} -F {} -B -s".format(bed, bed)
            else:
               flags = ""
            image = hgb_run(name_input, ref_id, range, coverage, flags, split, y, callet)
            stroke_color = st.sidebar.beta_color_picker("Stroke color hex: ")
            drawing_mode = st.sidebar.selectbox(
                "Drawing tool:", ("freedraw", "line", "rect", "circle", "transform"), 4
            )
            canvas_result = st_canvas(
                fill_color="rgba(255, 165, 0, 0.3)",  # Fixed fill color with some opacity
                stroke_width=3, #stroke_width,
                stroke_color=stroke_color,
                background_color="", #if bg_image else bg_color,
                background_image=image, 
                update_streamlit=False, #realtime_update or update_button,
                height=image.height,
                width=image.width,
                drawing_mode=drawing_mode,
                key="canvas",
            )             

        fields = ['seqid', 'start', 'end', 'source', 'featuretype', 'strand', 'attributes']
        allFoo =  list(db.region(region=(ref_id, range[0], range[1]), completely_within=False))
        df = pd.DataFrame([{fn: getattr(f, fn) for fn in fields} for f in allFoo], columns=fields)
        st.dataframe(df)
        #df = pd.DataFrame([vars(f) for f in list(db.region(region=(ref_id, range[0], range[1]), completely_within=False))])
        #df = pd.DataFrame(list(db.region(region=(ref_id, range[0], range[1]), completely_within=False)))
        #print(df)

    st.markdown(
        f"""
<style>
    .reportview-container .main .block-container{{
        max-width: 1280px;
    }}
</style>
""",
        unsafe_allow_html=True,
    )

