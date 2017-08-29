//C-  -*- C++ -*-
//C- -------------------------------------------------------------------
//C- DjVuLibre-3.5
//C- Copyright (c) 2002  Leon Bottou and Yann Le Cun.
//C- Copyright (c) 2001  AT&T
//C-
//C- This software is subject to, and may be distributed under, the
//C- GNU General Public License, either Version 2 of the license,
//C- or (at your option) any later version. The license should have
//C- accompanied the software or you may obtain a copy of the license
//C- from the Free Software Foundation at http://www.fsf.org .
//C-
//C- This program is distributed in the hope that it will be useful,
//C- but WITHOUT ANY WARRANTY; without even the implied warranty of
//C- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//C- GNU General Public License for more details.
//C- 
//C- DjVuLibre-3.5 is derived from the DjVu(r) Reference Library from
//C- Lizardtech Software.  Lizardtech Software has authorized us to
//C- replace the original DjVu(r) Reference Library notice by the following
//C- text (see doc/lizard2002.djvu and doc/lizardtech2007.djvu):
//C-
//C-  ------------------------------------------------------------------
//C- | DjVu (r) Reference Library (v. 3.5)
//C- | Copyright (c) 1999-2001 LizardTech, Inc. All Rights Reserved.
//C- | The DjVu Reference Library is protected by U.S. Pat. No.
//C- | 6,058,214 and patents pending.
//C- |
//C- | This software is subject to, and may be distributed under, the
//C- | GNU General Public License, either Version 2 of the license,
//C- | or (at your option) any later version. The license should have
//C- | accompanied the software or you may obtain a copy of the license
//C- | from the Free Software Foundation at http://www.fsf.org .
//C- |
//C- | The computer code originally released by LizardTech under this
//C- | license and unmodified by other parties is deemed "the LIZARDTECH
//C- | ORIGINAL CODE."  Subject to any third party intellectual property
//C- | claims, LizardTech grants recipient a worldwide, royalty-free, 
//C- | non-exclusive license to make, use, sell, or otherwise dispose of 
//C- | the LIZARDTECH ORIGINAL CODE or of programs derived from the 
//C- | LIZARDTECH ORIGINAL CODE in compliance with the terms of the GNU 
//C- | General Public License.   This grant only confers the right to 
//C- | infringe patent claims underlying the LIZARDTECH ORIGINAL CODE to 
//C- | the extent such infringement is reasonably necessary to enable 
//C- | recipient to make, have made, practice, sell, or otherwise dispose 
//C- | of the LIZARDTECH ORIGINAL CODE (or portions thereof) and not to 
//C- | any greater extent that may be necessary to utilize further 
//C- | modifications or combinations.
//C- |
//C- | The LIZARDTECH ORIGINAL CODE is provided "AS IS" WITHOUT WARRANTY
//C- | OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
//C- | TO ANY WARRANTY OF NON-INFRINGEMENT, OR ANY IMPLIED WARRANTY OF
//C- | MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
//C- +------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#if NEED_GNUG_PRAGMAS
# pragma implementation
#endif

/** @name cjb2

    {\bf Synopsis}
    \begin{verbatim}
        cjb2 [options] <input-pbm-or-tiff>  <output-djvu>
    \end{verbatim}

    {\bf Description}

    File #"cjb2.cpp"# demonstrates a simple encoder for Bilevel DjVu Images.
    It is able to perform lossless encoding and limited lossy encoding.  Lots
    of lossy encoding refinements are missing from this simple implementation.
    Comments in the code suggest a few improvements.

    Options are:
    \begin{description}
    \item[-dpi xxx]     Specify image resolution (default 300).
    \item[-lossless]    Lossless compression (same as -losslevel 0, default).
    \item[-clean]       Quasi-lossless compression (same as -losslevel 1).
    \item[-lossy]       Lossy compression (same as -losslevel 100).
    \item[-losslevel n] Set loss level (0 to 200)
    \item[-verbose]     Display additional messages.
    \end{description}
    Encoding is lossless unless one or several lossy options are selected.
    The #dpi# argument mostly affects the cleaning thresholds.

    {\bf Bugs}

    This is not the full-fledged multipage DjVu compressor, but merely a free
    tool provided with the DjVu Reference Library as a demonstrative example.

    @memo
    Simple JB2 encoder.
    @author
    L\'eon Bottou <leonb@research.att.com>\\
    Paul Howard <pgh@research.att.com>\\
    Pascal Vincent <vincentp@iro.umontreal.ca>\\
    Ilya Mezhirov <ilya@mezhirov.mccme.ru>
*/
//@{
//@}


#include "DjVuGlobal.h"
#include "GException.h"
#include "GSmartPointer.h"
#include "GContainer.h"
#include "ByteStream.h"
#include "IFFByteStream.h"
#include "GRect.h"
#include "GBitmap.h"
#include "JB2Image.h"
#include "DjVuInfo.h"
#include "GOS.h"
#include "GURL.h"
#include "DjVuMessage.h"
#include "jb2tune.h"
#include "common.h"
#if HAVE_TIFF
#include <tiffio.h>
#endif

#include <stdio.h>
#include <string>
#include <vector>

// --------------------------------------------------
// UTILITIES
// --------------------------------------------------

#ifdef MIN
#undef MIN
#endif

inline int 
MIN(int a, int b) 
{ 
  return ( a<b ?a :b); 
}

#ifdef MAX
#undef MAX
#endif

inline int 
MAX(int a, int b) 
{ 
  return ( a>b ?a :b); 
}






// --------------------------------------------------
// CONNECTED COMPONENT ANALYSIS AND CLEANING
// --------------------------------------------------

// -- A run of black pixels
struct Run    
{ 
  int y;         // vertical coordinate
  short x1;      // first horizontal coordinate
  short x2;      // last horizontal coordinate
  int ccid;      // component id
};

// -- A component descriptor
struct CC    
{
  GRect bb;      // bounding box
  int npix;      // number of black pixels
  int nrun;      // number of runs
  int frun;      // first run in cc ordered array of runs
};


// -- An image composed of runs
class CCImage 
{
public:
  int height;            // Height of the image in pixels
  int width;             // Width of the image in pixels
  GTArray<Run> runs;     // array of runs
  GTArray<CC>  ccs;      // Array of component descriptors
  int nregularccs;       // Number of regular ccs (set by merge_and_split_ccs)
  int largesize;         // CCs larger than that are special
  int smallsize;         // CCs smaller than that are special 
  int tinysize;          // CCs smaller than that may be removed 
  CCImage();
  void init(int width, int height, int dpi);
  void add_single_run(int y, int x1, int x2, int ccid=0);
  void add_bitmap_runs(const GBitmap &bm, int offx=0, int offy=0, int ccid=0);
  GP<GBitmap> get_bitmap_for_cc(int ccid) const;
  void append_to_jb2image(JB2Image* jimg) const;
  void make_ccids_by_analysis();
  void make_ccs_from_ccids();
  void erase_tiny_ccs();
  void merge_and_split_ccs();
  void sort_in_reading_order(); 
};


// -- Compares runs
static inline bool
operator <= (const Run &a, const Run &b)
{
  return (a.y<b.y) || (a.y==b.y && a.x1<=b.x1);
}


// -- Constructs CCImage and provide defaults
CCImage::CCImage()
  : height(0), width(0), nregularccs(0)
{
}

void
CCImage::init(int w, int h, int dpi)
{
  runs.empty();
  ccs.empty();
  height = h;
  width = w;
  nregularccs = 0;
  dpi = MAX(200, MIN(900, dpi));
  largesize = MIN( 500, MAX(64, dpi));
  smallsize = MAX(2, dpi/150);
  tinysize = MAX(0, dpi*dpi/20000 - 1);
}


// -- Adds a run to the CCImage
inline void 
CCImage::add_single_run(int y, int x1, int x2, int ccid)
{
  int index = runs.hbound();
  runs.touch(++index);
  Run& run = runs[index];
  run.y = y;
  run.x1 = x1;
  run.x2 = x2;
  run.ccid = ccid;
}


// -- Adds runs extracted from a bitmap
void 
CCImage::add_bitmap_runs(const GBitmap &bm, int offx, int offy, int ccid)
{
  // Iterate over rows
  for (unsigned int y=0; y<bm.rows(); y++)
    {
      const unsigned char *row = bm[y];
      int w = bm.columns();
      int x = 0;
      // Iterate over runs
      while (x < w)
        {
          while (x < w  && !row[x]) x++;
          if (x < w)
            {
              int x1 = x;
              while (x < w && row[x]) x++;
              add_single_run(offy+y, offx+x1, offx+x-1, ccid);
            }
        }
    }
}


// -- Performs connected component analysis
void
CCImage::make_ccids_by_analysis()
{
  // Sort runs
  runs.sort();
  // Single Pass Connected Component Analysis (with unodes)
  int n;
  int p=0;
  GTArray<int> umap;
  for (n=0; n<=runs.hbound(); n++)
    {
      int y = runs[n].y;
      int x1 = runs[n].x1 - 1;
      int x2 = runs[n].x2 + 1;
      int id = (umap.hbound() + 1);
      // iterate over previous line runs
      for(;runs[p].y < y-1;p++);
      for(;(runs[p].y < y) && (runs[p].x1 <= x2);p++ )
        {
          if ( runs[p].x2 >= x1 )
            {
              // previous run touches current run
              int oid = runs[p].ccid;
              while (umap[oid] < oid)
                oid = umap[oid];
              if ((int)id > umap.hbound()) {
                id = oid;
              } else if (id < oid) {
                umap[oid] = id;
              } else {
                umap[id] = oid;
                id = oid;
              }
              // freshen previous run id
              runs[p].ccid = id;
              // stop if previous run goes past current run
              if (runs[p].x2 >= x2)
                break;
            }
        }
      // create new entry in umap
      runs[n].ccid = id;
      if (id > umap.hbound())
        {
          umap.touch(id);
          umap[id] = id;
        }
    }
  // Update umap and ccid
  for (n=0; n<=runs.hbound(); n++)
    {
      Run &run = runs[n];
      int ccid = run.ccid;
      while (umap[ccid] < ccid)
      {
        ccid = umap[ccid];
      }
      umap[run.ccid] = ccid;
      run.ccid = ccid;
    }
}


// -- Constructs the ``ccs'' array from run's ccids.
void
CCImage::make_ccs_from_ccids()
{
  int n;
  Run *pruns = runs;
  // Find maximal ccid
  int maxccid = nregularccs-1;
  for (n=0; n<=runs.hbound(); n++)
    if (pruns[n].ccid > maxccid)
      maxccid = runs[n].ccid;
  // Renumber ccs 
  GTArray<int> armap(0,maxccid);
  int *rmap = armap;
  for (n=0; n<=maxccid; n++)
    armap[n] = -1;
  for (n=0; n<=runs.hbound(); n++)
    if (pruns[n].ccid >= 0)
      rmap[ pruns[n].ccid ] = 1;
  int nid = 0;
  for (n=0; n<=maxccid; n++)
    if (rmap[n] > 0)
      rmap[n] = nid++;
  
  // Adjust nregularccs (since ccs are renumbered)
  while (nregularccs>0 && rmap[nregularccs-1]<0)
    nregularccs -= 1;
  if (nregularccs>0)
    nregularccs = 1 + rmap[nregularccs-1];
  
  // Prepare cc descriptors
  ccs.resize(0,nid-1);
  for (n=0; n<nid; n++)
    ccs[n].nrun = 0;
  
  // Relabel runs
  for (n=0; n<=runs.hbound(); n++)
    {
      Run &run = pruns[n];
      if (run.ccid < 0) continue;  // runs with negative ccids are destroyed
      int oldccid = run.ccid;
      int newccid = rmap[oldccid];
      CC &cc = ccs[newccid];
      run.ccid = newccid;
      cc.nrun += 1;
    }
  
  // Compute positions for runs of cc
  int frun = 0;
  for (n=0; n<nid; n++) 
    {
      ccs[n].frun = rmap[n] = frun;
      frun += ccs[n].nrun;
    }

  // Copy runs
  GTArray<Run> rtmp;
  rtmp.steal(runs);
  Run *ptmp = rtmp;
  runs.resize(0,frun-1);
  pruns = runs;
  for (n=0; n<=rtmp.hbound(); n++)
    {
      int id = ptmp[n].ccid;
      if (id < 0) continue;
      int pos = rmap[id]++;
      pruns[pos] = ptmp[n];
    }

  // Finalize ccs
  for (n=0; n<nid; n++)
    {
      CC &cc = ccs[n];
      int npix = 0;
      runs.sort(cc.frun, cc.frun+cc.nrun-1);
      Run *run = &runs[cc.frun];
      int xmin = run->x1;
      int xmax = run->x2;
      int ymin = run->y;
      int ymax = run->y;
      for (int i=0; i<cc.nrun; i++, run++)
        {
          if (run->x1 < xmin)  xmin = run->x1;
          if (run->x2 > xmax)  xmax = run->x2;
          if (run->y  < ymin)  ymin = run->y;
          if (run->y  > ymax)  ymax = run->y;
          npix += run->x2 - run->x1 + 1;
        }
      cc.npix = npix;
      cc.bb.xmin = xmin;
      cc.bb.ymin = ymin;
      cc.bb.xmax = xmax + 1;
      cc.bb.ymax = ymax + 1;
    }
}


// Removes ccs which are too small.
void
CCImage::erase_tiny_ccs()
{
  // ISSUE: HALFTONE DETECTION
  // We should not remove tiny ccs if they are part of a halftone pattern...
  for (int i=0; i<ccs.size(); i++)
    {
      CC& cc = ccs[i];
      if (cc.npix <= tinysize)
        {
          // Mark cc to be erased
          Run *r = &runs[cc.frun];
          int nr = cc.nrun;
          cc.nrun = 0;
          cc.npix = 0;
          while (--nr >= 0)
            (r++)->ccid = -1;
        }
    }
}
 

// -- Merges small ccs and split large ccs
void
CCImage::merge_and_split_ccs()
{
  int ncc = ccs.size();
  int nruns = runs.size();
  int splitsize = largesize;
  if (ncc <= 0) return;
  // Grid of special components
  int gridwidth = (width+splitsize-1)/splitsize;
  nregularccs = ncc;
  // Set the correct ccids for the runs
  for (int ccid=0; ccid<ncc; ccid++)
    {
      CC* cc = &ccs[ccid];
      if (cc->nrun <= 0) continue;
      int ccheight = cc->bb.height();
      int ccwidth = cc->bb.width();
      if (ccheight<=smallsize && ccwidth<=smallsize)
        {
          int gridi = (cc->bb.ymin+cc->bb.ymax)/splitsize/2;
          int gridj = (cc->bb.xmin+cc->bb.xmax)/splitsize/2;
          int newccid = ncc + gridi*gridwidth + gridj;
          for(int runid=cc->frun; runid<cc->frun+cc->nrun; runid++)
            runs[runid].ccid = newccid;
        }
      else if (ccheight>=largesize || ccwidth>=largesize)
        {
          for(int runid=cc->frun; runid<cc->frun+cc->nrun; runid++)
            {
              Run& r = runs[runid];
              int y = r.y;
              int x_start = r.x1;
              int x_end = r.x2;
              int gridi = y/splitsize;
              int gridj_start = x_start/splitsize;
              int gridj_end = x_end/splitsize;
              int gridj_span = gridj_end-gridj_start;
              int newccid = ncc + gridi*gridwidth + gridj_start;
              if (! gridj_span)
                {
                  r.ccid = newccid;
                }
              else // gridj_span>0
                {
                  // truncate the current run 
                  r.ccid = newccid++;
                  int x = (gridj_start+1)*splitsize;
                  r.x2 = x-1;
                  runs.touch(nruns+gridj_span-1);
                  // append additional runs to the runs array
                  for(int gridj=gridj_start+1; gridj<gridj_end; gridj++)
                    {
                      Run& newrun = runs[nruns++];
                      newrun.y = y;
                      newrun.x1 = x;
                      x += splitsize;
                      newrun.x2 = x-1;
                      newrun.ccid = newccid++;
                    }
                  // append last run to the run array
                  Run& newrun = runs[nruns++];
                  newrun.y = y;
                  newrun.x1 = x;
                  newrun.x2 = x_end;
                  newrun.ccid = newccid++;                      
                }
            }
        }
    }
  // Recompute cc descriptors
  make_ccs_from_ccids();
}


// -- Helps sorting cc
static int 
top_edges_descending (const void *pa, const void *pb)
{
  if (((CC*) pa)->bb.ymax != ((CC*) pb)->bb.ymax)
    return (((CC*) pb)->bb.ymax - ((CC*) pa)->bb.ymax);
  if (((CC*) pa)->bb.xmin != ((CC*) pb)->bb.xmin)
    return (((CC*) pa)->bb.xmin - ((CC*) pb)->bb.xmin);
  return (((CC*) pa)->frun - ((CC*) pb)->frun);
}


// -- Helps sorting cc
static int 
left_edges_ascending (const void *pa, const void *pb)
{
  if (((CC*) pa)->bb.xmin != ((CC*) pb)->bb.xmin)
    return (((CC*) pa)->bb.xmin - ((CC*) pb)->bb.xmin);
  if (((CC*) pb)->bb.ymax != ((CC*) pa)->bb.ymax)
    return (((CC*) pb)->bb.ymax - ((CC*) pa)->bb.ymax);
  return (((CC*) pa)->frun - ((CC*) pb)->frun);
}


// -- Helps sorting cc
static int 
integer_ascending (const void *pa, const void *pb)
{
  return ( *(int*)pb - *(int*)pa );
}


// -- Sort ccs in approximate reading order
void 
CCImage::sort_in_reading_order()
{
  if (nregularccs<2) return;
  CC *ccarray = new CC[nregularccs];
  // Copy existing ccarray (but segregate special ccs)
  int ccid;
  for(ccid=0; ccid<nregularccs; ccid++)
    ccarray[ccid] = ccs[ccid];
  // Sort the ccarray list into top-to-bottom order.
  qsort (ccarray, nregularccs, sizeof(CC), top_edges_descending);
  // Subdivide the ccarray list roughly into text lines [LYB]
  // - Determine maximal top deviation
  int maxtopchange = width / 40;
  if (maxtopchange < 32) 
    maxtopchange = 32;
  // - Loop until processing all ccs
  int ccno = 0;
  int *bottoms = new int[nregularccs];
  while (ccno < nregularccs)
    {
      // - Gather first line approximation
      int nccno;
      int sublist_top = ccarray[ccno].bb.ymax-1;
      int sublist_bottom = ccarray[ccno].bb.ymin;
      for (nccno=ccno; nccno < nregularccs; nccno++)
        {
          if (ccarray[nccno].bb.ymax-1 < sublist_bottom) break;
          if (ccarray[nccno].bb.ymax-1 < sublist_top - maxtopchange) break;
          int bottom = ccarray[nccno].bb.ymin;
          bottoms[nccno-ccno] = bottom;
          if (bottom < sublist_bottom)
            sublist_bottom = bottom;
        }
      // - If more than one candidate cc for the line
      if (nccno > ccno + 1)
        {
          // - Compute median bottom
          qsort(bottoms, nccno-ccno, sizeof(int), integer_ascending);
          int bottom = bottoms[ (nccno-ccno-1)/2 ];
          // - Compose final line
          for (nccno=ccno; nccno < nregularccs; nccno++)
            if (ccarray[nccno].bb.ymax-1 < bottom)
              break;
          // - Sort final line
          qsort (ccarray+ccno, nccno-ccno, sizeof(CC), left_edges_ascending);
        }
      // - Next line
      ccno = nccno;
    }
  // Copy ccarray back and renumber the runs
  for(ccid=0; ccid<nregularccs; ccid++)
    {
      CC& cc = ccarray[ccid];
      ccs[ccid] = cc;
      for(int r=cc.frun; r<cc.frun+cc.nrun; r++)
        runs[r].ccid = ccid;
    }
  // Free memory
  delete [] bottoms;
  delete[] ccarray;
}


// -- Creates a bitmap for a particular component
GP<GBitmap>   
CCImage::get_bitmap_for_cc(const int ccid) const
{
  const CC &cc = ccs[ccid];
  const GRect &bb = cc.bb;
  GP<GBitmap> bits = GBitmap::create(bb.height(), bb.width());
  const Run *prun = & runs[(int)cc.frun];
  for (int i=0; i<cc.nrun; i++,prun++)
    {
      if (prun->y<bb.ymin || prun->y>=bb.ymax)
        G_THROW("Internal error (y bounds)");
      if (prun->x1<bb.xmin || prun->x2>=bb.xmax)
        G_THROW("Internal error (x bounds)");
      unsigned char *row = (*bits)[prun->y - bb.ymin];
      for (int x=prun->x1; x<=prun->x2; x++)
        row[x - bb.xmin] = 1;
    }
  return bits;
}


// -- Append the remaining components to a JB2Image
void
CCImage::append_to_jb2image(JB2Image *jimg) const
{
  if (runs.hbound() < 0) return;
  if (ccs.hbound() < 0)
    G_THROW("Must first perform a cc analysis");
  // Iterate over CCs
  for (int ccid=0; ccid<=ccs.hbound(); ccid++)
    {
      JB2Shape shape;
      JB2Blit  blit;
      shape.parent = -1;
      shape.bits = get_bitmap_for_cc(ccid);
      shape.userdata = 0;
      if (ccid >= nregularccs)
        shape.userdata |= JB2SHAPE_SPECIAL;
      blit.shapeno = jimg->add_shape(shape);
      blit.left = ccs[ccid].bb.xmin;
      blit.bottom = ccs[ccid].bb.ymin;
      jimg->add_blit(blit);
      shape.bits->compress();
    }
  // Return
  return;
}




// --------------------------------------------------
// COMPLETE COMPRESSION ROUTINE
// --------------------------------------------------

struct cjb2opts {
  int  dpi;
  int  forcedpi;
  int  losslevel;
  bool verbose;
  GURL dict;
  std::string output_dict;
  int pages_per_dict;
};

#if HAVE_TIFF

static int 
is_tiff(ByteStream *ref) 
{
  char magic[2];
  magic[0] = magic[1] = 0;
  ref->readall((void*)magic, sizeof(magic));
  ref->seek(0);
  if(magic[0] == 'I' && magic[1] == 'I')
    return 1;
  if(magic[0] == 'M' && magic[1] == 'M')
    return 1;
  return 0;
}

static tsize_t readproc(thandle_t h, tdata_t p, tsize_t s) { 
  ByteStream *bs = (ByteStream*)h;
  return (tsize_t) bs->readall((void*)p, (size_t)s);
}

static tsize_t writeproc(thandle_t, tdata_t, tsize_t) { 
  return -1; 
}

static toff_t seekproc(thandle_t h, toff_t offset, int mode) { 
  ByteStream *bs = (ByteStream*)h;
  bs->seek((long)offset, mode);
  return (toff_t)bs->tell(); 
}

static int closeproc(thandle_t) { 
  return 0; 
}

static toff_t sizeproc(thandle_t h) { 
  ByteStream *bs = (ByteStream*)h;
  return (toff_t) bs->size(); 
}

static int mapproc(thandle_t, tdata_t*, toff_t*) { 
  return -1; 
}

static void unmapproc(thandle_t, tdata_t, toff_t) { 
}

static int get_tiff_page_count(ByteStream *bs) {
	bs->seek(0);
	TIFF *tiff = TIFFClientOpen("libtiff", "rm", (thandle_t)bs,
		readproc, writeproc, seekproc,
		closeproc, sizeproc,
		mapproc, unmapproc);
	int dircount = 0;
	if (tiff) {
		do {
			dircount++;
		} while (TIFFReadDirectory(tiff));
		TIFFClose(tiff);
	}
	return dircount;
}

static TIFF* open_tiff(ByteStream *bs) {
	bs->seek(0);
	return TIFFClientOpen("libtiff", "rm", (thandle_t)bs,
		readproc, writeproc, seekproc,
		closeproc, sizeproc,
		mapproc, unmapproc);
}

static void
read_tiff(CCImage &rimg, TIFF *tiff, cjb2opts &opts)
{
  // bitonal
  uint16 bps = 0, spp = 0;
  TIFFGetFieldDefaulted(tiff, TIFFTAG_BITSPERSAMPLE, &bps);
  TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLESPERPIXEL, &spp);
  if (bps != 1 || spp != 1)
    G_THROW("Tiff image is not bitonal");
  // photometric
  uint16 photo = PHOTOMETRIC_MINISWHITE;
  TIFFGetFieldDefaulted(tiff, TIFFTAG_PHOTOMETRIC, &photo);
  // image size
  uint32 w, h;
  if (!TIFFGetFieldDefaulted(tiff, TIFFTAG_IMAGEWIDTH, &w) ||
      !TIFFGetFieldDefaulted(tiff, TIFFTAG_IMAGELENGTH, &h) )
    G_THROW("Tiff image size is not defined");
  // resolution
  float xres, yres;
  if (TIFFGetFieldDefaulted(tiff, TIFFTAG_XRESOLUTION, &xres) &&
      TIFFGetFieldDefaulted(tiff, TIFFTAG_YRESOLUTION, &yres) ) 
    {
      if (xres != yres)
        DjVuPrintErrorUTF8( "cjb2: X- and Y-resolution do not match\n");
      if (! opts.forcedpi)
        opts.dpi = (int) (xres + yres) / 2;
    }
  // init rimg
  rimg.init(w, h, opts.dpi);
  // allocate scanline
  tsize_t scanlinesize = TIFFScanlineSize(tiff);
  scanlinesize = MAX(scanlinesize,1);
  unsigned char *scanline = 0;
  GPBuffer<unsigned char> gscanline(scanline, scanlinesize);
  // iterate on rows
  for (int y=0; y<(int)h; y++)
    {
      int yy = h - y - 1;
      if (TIFFReadScanline(tiff, (tdata_t)scanline, y) < 0)
        G_THROW("Tiff file is corrupted (TIFFReadScanline)");
      if (photo != PHOTOMETRIC_MINISWHITE)
        for (int i=0; i<(int)scanlinesize; i++)
          scanline[i] ^= 0xff;
      int lastx=0, off=0;
      unsigned char mask=0, c=0, b=0;
      for (int x=0; x<(int)w; x++)
        {
          if (! mask) 
            {
              b = scanline[off++];
              while (b==c && x+8<(int)w )
                {
                  x = x + 8;  // speedup
                  b = scanline[off++];
                }
              mask = 0x80;
            }
          if ( (b ^ c) & mask ) 
            {
              c ^= 0xff;
              if (c)
                lastx = x;
              else
                rimg.add_single_run(yy, lastx, x-1);
            }
          mask >>= 1;
        }
      if (c)
        rimg.add_single_run(yy, lastx, w-1);
    }
}

#endif // HAVE_TIFF

// Load dictionary recursively
GP<JB2Dict> loadDictionary(const GURL& fileName, int recursive = 16) {
	GP<JB2Dict> dict = JB2Dict::create();
	GP<ByteStream> ibs = ByteStream::create(fileName, "rb");
	GP<IFFByteStream> iiff = IFFByteStream::create(ibs);
	GUTF8String chkid;
	iiff->get_chunk(chkid);
	if (chkid != "FORM:DJVI") {
		DjVuFormatErrorUTF8("The file '%s' is not FORM:DJVI file\n", (const char*)fileName);
	} else {
		for (;;) {
			chkid.empty();
			if (iiff->get_chunk(chkid) == 0 || chkid.length() == 0) {
				DjVuFormatErrorUTF8("The file '%s' does not contain Djbz chunk\n", (const char*)fileName);
				break;
			}
			if (chkid == "INCL") {
				if (recursive <= 0) {
					DjVuFormatErrorUTF8("Maximal recursive count reached; give up\n");
					break;
				}

				GP<ByteStream> ibs2 = iiff->get_bytestream();
				const int M = 1024;
				int m = ibs2->size();
				if (m <= 0 || m >= M) {
					DjVuFormatErrorUTF8("INCL size too big or too small\n");
					break;
				}
				char s[M];
				ibs2->readall(s, m);
				s[m] = '\0';

				GP<JB2Dict> dict2 = loadDictionary(GURL::Filename::UTF8(s), recursive - 1);
				if (dict2) {
					dict->set_inherited_dict(dict2);
				} else {
					break;
				}
			} else if (chkid == "Djbz") {
				dict->decode(iiff->get_bytestream());
				return dict;
			}
			iiff->close_chunk();
		}
	}
	return 0;
}

class cjb2 {
public:
	cjb2(const std::vector<GURL> &inputlist, const std::string &outputname_, cjb2opts &opts_);
private:
	int pageno;
	int dictno;
	bool isMultipage;
	bool isMultiDict;
	cjb2opts &opts;
	GP<JB2Dict> shared_dict;
	GP<JB2Image> jimg; // image which contains data of multiple pages
	std::vector<int> jimg_width;
	std::vector<int> jimg_height;
	std::vector<int> jimg_shapes;
	std::vector<int> jimg_blits;
	std::string outputname;

	void cjb2_output();
	GURL get_output_name();
	GURL get_output_dict_name();
	void print_tune_result();
};

GURL cjb2::get_output_name() {
	GURL urlout;
	if (isMultipage) {
		int M = outputname.size() + 1024;
		char *s = new char[M];
#ifdef WIN32
		_snprintf(s, M, outputname.c_str(), pageno);
#else
		snprintf(s, M, outputname.c_str(), pageno);
#endif
		s[M - 1] = '\0';
		urlout = GURL::Filename::UTF8(s);
		delete[] s;
	} else {
		urlout = GURL::Filename::UTF8(outputname.c_str());
	}
	return urlout;
}

GURL cjb2::get_output_dict_name() {
	GURL urlout;
	if (isMultiDict) {
		int M = opts.output_dict.size() + 1024;
		char *s = new char[M];
#ifdef WIN32
		_snprintf(s, M, outputname.c_str(), dictno);
#else
		snprintf(s, M, outputname.c_str(), dictno);
#endif
		s[M - 1] = '\0';
		urlout = GURL::Filename::UTF8(s);
		delete[] s;
	} else {
		urlout = GURL::Filename::UTF8(opts.output_dict.c_str());
	}
	return urlout;
}

void cjb2::print_tune_result() {
	if (opts.verbose)
	{
		int nshape = 0, nrefine = 0;
		for (int i = jimg->get_inherited_shape_count(), m = jimg->get_shape_count(); i<m; i++) {
			if (!jimg->get_shape(i).bits) continue;
			if (jimg->get_shape(i).parent >= 0) nrefine++;
			nshape++;
		}
		DjVuFormatErrorUTF8("%s\t%d\t%d", ERR_MSG("cjb2.shapes"),
			nshape, nrefine);
	}
}

void cjb2::cjb2_output() {
	if (!jimg || jimg_width.empty()) return; // do nothing

	// Pattern matching
	if (opts.losslevel>1)
		tune_jb2image_lossy(jimg, opts.dpi, opts.losslevel);
	else
		tune_jb2image_lossless(jimg);

	print_tune_result();

	if (opts.pages_per_dict <= 1) { // does not generate dictionary at all
		jimg->set_dimension(jimg_width[0], jimg_height[0]);

		// get the output file name
		pageno++;
		GURL urlout = get_output_name();

		// Code
		GP<ByteStream> obs = ByteStream::create(urlout, "wb");
		GP<IFFByteStream> giff = IFFByteStream::create(obs);
		IFFByteStream &iff = *giff;
		// -- main composite chunk
		iff.put_chunk("FORM:DJVU", 1);
		// -- ``INFO'' chunk
		GP<DjVuInfo> ginfo = DjVuInfo::create();
		DjVuInfo &info = *ginfo;
		info.height = jimg_height[0];
		info.width = jimg_width[0];
		info.dpi = opts.dpi;
		iff.put_chunk("INFO");
		info.encode(*iff.get_bytestream());
		iff.close_chunk();
		if (jimg->get_blit_count() > 0) {
			// -- ``INCL'' chunk
			if (shared_dict) {
				iff.put_chunk("INCL");
				GUTF8String fname = opts.dict.fname();
				const char* s = fname;
				iff.get_bytestream()->writall(s, strlen(s));
				iff.close_chunk();
			}
			// -- ``Sjbz'' chunk
			iff.put_chunk("Sjbz");
			jimg->encode(iff.get_bytestream());
			iff.close_chunk();
		}
		// -- terminate main composite chunk
		iff.close_chunk();
		// Finished!
	} else {
		int shapes0 = jimg->get_inherited_shape_count();
		int shapes = jimg->get_shape_count();
		int shapes1 = shapes - shapes0;
		int shapes0new = 0;

		std::vector<int> shapeFreq, newIndex;

		GP<JB2Dict> newDict = JB2Dict::create();
		if (shared_dict) newDict->set_inherited_dict(shared_dict);

		if (shapes1 > 0) {
			// count the frequency
			shapeFreq.resize(shapes1, 0);
			newIndex.resize(shapes1, -1);
			for (int i = 0, m = jimg->get_blit_count(); i < m; i++) {
				int shapeno = jimg->get_blit(i)->shapeno;
				if (shapeno >= shapes0 && jimg->get_shape(shapeno).bits) shapeFreq[shapeno - shapes0]++;
			}
			for (int i = shapes0; i < shapes; i++) {
				JB2Shape &shape = jimg->get_shape(i);
				if (shape.bits) {
					int shapeno = shape.parent;
					if (shapeno >= shapes0 && jimg->get_shape(shapeno).bits) shapeFreq[shapeno - shapes0] += 2;
				}
			}

			// create the new dictionary, calculate the new index
			for (int i = 0; i < shapes1; i++) {
				if (shapeFreq[i] > 1) {
					JB2Shape shape = jimg->get_shape(shapes0 + i);
					if (shape.parent >= shapes0) { // relocate
						if ((shape.parent = newIndex[shape.parent - shapes0]) < 0) {
							G_THROW("BUG: the new index should be valid");
						}
					}
					newIndex[i] = newDict->add_shape(shape);
					shapes0new++;
				}
			}
		}

		// save new dictionary
		dictno++;
		GURL dictout = get_output_dict_name();
		{
			GP<ByteStream> obs = ByteStream::create(dictout, "wb");
			GP<IFFByteStream> giff = IFFByteStream::create(obs);
			IFFByteStream &iff = *giff;
			// -- main composite chunk
			iff.put_chunk("FORM:DJVI", 1);
			// -- ``INCL'' chunk
			if (shared_dict) {
				iff.put_chunk("INCL");
				GUTF8String fname = opts.dict.fname();
				const char* s = fname;
				iff.get_bytestream()->writall(s, strlen(s));
				iff.close_chunk();
			}
			// -- ``Djbz'' chunk
			iff.put_chunk("Djbz");
			newDict->encode(iff.get_bytestream());
			iff.close_chunk();
			// -- terminate main composite chunk
			iff.close_chunk();
		}

		// for each page
		int shapeno = shapes0;
		int blitno = 0;
		for (int i = 0, m = jimg_width.size(); i < m; i++) {
			// create new JB2Image
			GP<JB2Image> jimg2 = JB2Image::create();
			jimg2->set_inherited_dict(newDict);
			jimg2->set_dimension(jimg_width[i], jimg_height[i]);

			// calculate the new index for remaining shapes
			for (int m2 = jimg_shapes[i]; shapeno < m2; shapeno++) {
				if (shapeFreq[shapeno - shapes0] == 1) {
					JB2Shape shape = jimg->get_shape(shapeno);
					if (shape.parent >= shapes0) { // relocate
						if ((shape.parent = newIndex[shape.parent - shapes0]) < 0) {
							G_THROW("BUG: the new index should be valid");
						}
					}
					int idx = jimg2->add_shape(shape);
					if (newIndex[shapeno - shapes0] >= 0) {
						G_THROW("BUG: the new index should be -1");
					}
					newIndex[shapeno - shapes0] = idx;
				}
			}

			// put the blits to new JB2Image
			for (int m2 = jimg_blits[i]; blitno < m2; blitno++) {
				JB2Blit blit = *jimg->get_blit(blitno);
				if ((int)blit.shapeno >= shapes0) { // relocate
					blit.shapeno = newIndex[blit.shapeno - shapes0];
					if (blit.shapeno < 0) continue;
				}
				jimg2->add_blit(blit);
			}

			// save new file
			pageno++;
			GURL urlout = get_output_name();
			{
				GP<ByteStream> obs = ByteStream::create(urlout, "wb");
				GP<IFFByteStream> giff = IFFByteStream::create(obs);
				IFFByteStream &iff = *giff;
				// -- main composite chunk
				iff.put_chunk("FORM:DJVU", 1);
				// -- ``INFO'' chunk
				GP<DjVuInfo> ginfo = DjVuInfo::create();
				DjVuInfo &info = *ginfo;
				info.height = jimg_height[i];
				info.width = jimg_width[i];
				info.dpi = opts.dpi;
				iff.put_chunk("INFO");
				info.encode(*iff.get_bytestream());
				iff.close_chunk();
				if (jimg2->get_blit_count() > 0) {
					// -- ``INCL'' chunk
					iff.put_chunk("INCL");
					GUTF8String fname = dictout.fname();
					const char* s = fname;
					iff.get_bytestream()->writall(s, strlen(s));
					iff.close_chunk();
					// -- ``Sjbz'' chunk
					iff.put_chunk("Sjbz");
					jimg2->encode(iff.get_bytestream());
					iff.close_chunk();
				}
				// -- terminate main composite chunk
				iff.close_chunk();
			}
		}
	}

	// erase processed data
	jimg = 0;
	jimg_width.clear();
	jimg_height.clear();
	jimg_shapes.clear();
	jimg_blits.clear();
}

cjb2::cjb2(const std::vector<GURL> &inputlist, const std::string &outputname_, cjb2opts &opts_)
	: pageno(0)
	, dictno(0)
	, isMultipage(false)
	, isMultiDict(false)
	, opts(opts_)
	, outputname(outputname_)
{
	// Load shared dictionary
	if (!opts.dict.is_empty()) {
		shared_dict = loadDictionary(opts.dict);
	}

	// get the number of pages
	std::vector<int> pagecounts;
	for (int i = 0, m = inputlist.size(); i < m; i++) {
		int pagecount = 0;
#if HAVE_TIFF
		GP<ByteStream> ibs = ByteStream::create(inputlist[i], "rb");
		if (is_tiff(ibs)) {
			pagecount = get_tiff_page_count(ibs);
		} else {
			// assume all non-TIFF files are of single page
			pagecount = 1;
		}
#else
		// assume all files are of single page
		pagecount = 1;
#endif
		pagecounts.push_back(pagecount);
		pageno += pagecount;
	}

	if (pageno <= 0) return; // do nothing

	if (opts.output_dict.empty()) {
		opts.pages_per_dict = 1;
	} else {
		if (opts.pages_per_dict <= 0) opts.pages_per_dict = 0x7FFFFFF0;
	}
	
	if (pageno > 1) {
		isMultipage = true;
		// rename the output file to contain '%d'
		if (outputname.find_last_of('%') == std::string::npos) {
			size_t lpe = outputname.find_last_of('.');
			if (lpe == std::string::npos) {
				outputname += ".%04d.djvu";
			} else {
				outputname = outputname.substr(0, lpe) + ".%04d" + outputname.substr(lpe);
			}
		}

		if (opts.pages_per_dict > 1 && pageno > opts.pages_per_dict) {
			// we need multiple dictionaries
			isMultiDict = true;
			if (opts.output_dict.find_last_of('%') == std::string::npos) {
				size_t lpe = opts.output_dict.find_last_of('.');
				if (lpe == std::string::npos) {
					opts.output_dict += ".%04d";
				} else {
					opts.output_dict = opts.output_dict.substr(0, lpe) + ".%04d" + opts.output_dict.substr(lpe);
				}
			}
		}
	}

	// reset counter, etc.
	pageno = 0;
	dictno = 0;

	// for each file
	for (int i = 0, m = inputlist.size(); i < m; i++) {
		int pagecount = pagecounts[i];

		GP<ByteStream> ibs = ByteStream::create(inputlist[i], "rb");

		// open TIFF file
#if HAVE_TIFF
		TIFF *tiff = 0;
		if (is_tiff(ibs)) {
			tiff = open_tiff(ibs);
		}
#endif

		// for each page in the file
		for (int p = 0; p < pagecount; p++) {
			CCImage rimg;

			// load the current page
#if HAVE_TIFF
			if (tiff) {
				if (p > 0) {
					if(!TIFFReadDirectory(tiff)) break;
				}
				read_tiff(rimg, tiff, opts);
			} else
#endif
			{
				GP<GBitmap> input = GBitmap::create(*ibs);
				rimg.init(input->columns(), input->rows(), opts.dpi);
				rimg.add_bitmap_runs(*input);
			}
			if (opts.verbose)
				DjVuFormatErrorUTF8("%s\t%d", ERR_MSG("cjb2.runs"),
				rimg.runs.size());

			// Component analysis
			rimg.make_ccids_by_analysis(); // obtain ccids
			rimg.make_ccs_from_ccids();    // compute cc descriptors
			if (opts.verbose)
				DjVuFormatErrorUTF8("%s\t%d", ERR_MSG("cjb2.ccs_before"),
				rimg.ccs.size());
			if (opts.losslevel > 0)
				rimg.erase_tiny_ccs();       // clean
			rimg.merge_and_split_ccs();    // reorganize weird ccs
			rimg.sort_in_reading_order();  // sort cc descriptors
			if (opts.verbose)
				DjVuFormatErrorUTF8("%s\t%d", ERR_MSG("cjb2.ccs_after"),
				rimg.ccs.size());

			// Append to jb2image (assuming tune_jb2image_* doesn't change the order of shapes and blits)
			if (!jimg) {
				jimg = JB2Image::create();
				if (shared_dict) jimg->set_inherited_dict(shared_dict);
			}
			rimg.append_to_jb2image(jimg);
			jimg_width.push_back(rimg.width);
			jimg_height.push_back(rimg.height);
			jimg_shapes.push_back(jimg->get_shape_count());
			jimg_blits.push_back(jimg->get_blit_count());

			// save memory
			rimg.runs.empty();
			rimg.ccs.empty();

			// Check if we need to output result
			if ((int)jimg_width.size() >= opts.pages_per_dict) cjb2_output();
		}

		// close TIFF file
#if HAVE_TIFF
		if (tiff) {
			TIFFClose(tiff);
		}
#endif
	}

	// Output remaining result
	cjb2_output();
}





// --------------------------------------------------
// MAIN
// --------------------------------------------------

void
usage()
{
  DjVuPrintErrorUTF8(
#ifdef DJVULIBRE_VERSION
         "CJB2 --- DjVuLibre-" DJVULIBRE_VERSION "\n"
#endif
         "Simple DjVuBitonal encoder\n\n"
         "Usage: cjb2 [options] (<input-pbm-or-tiff> | -f <filelist> ) ... <output-djvu>\n"
         "Options are:\n"
         " -v, -verbose    Display additional messages.\n"
         " -dpi <n>        Specify image resolution (default 300).\n"
         " -clean          Cleanup image by removing small flyspecks.\n"
         " -lossy          Lossy compression (implies -clean as well)\n"
         " -losslevel <n>  Loss factor (implies -lossy, default 100)\n"
         " -dict <file>    Set input dictionary file (experimental)\n"
		 " -output-dict <file>\n"
		 "                 Generate dictionary for multipage encoding (experimental)\n"
		 " -p <n>, -pages-per-dict <n>\n"
		 "                 Pages per dictionary (default 10, 0=all)\n"
         "Encoding is lossless unless a lossy options is selected.\n" );
  exit(10);
}

static void add_filelist(const GURL& filename, std::vector<GURL>& list) {
	GP<ByteStream> f = ByteStream::create(filename, "rb");

	std::string s;

	for (;;) {
		int c;
		unsigned char ch;
		if (f->readall(&ch, 1) == 1)
			c = ch;
		else
			c = -1;

		if (c == '\r') continue;
		else if (c < 0 || c == '\n') {
			size_t pos = s.find_first_not_of(" \t");
			if (pos != s.npos) s = s.substr(pos);
			else s.clear();
			pos = s.find_last_not_of(" \t");
			if (pos != s.npos) s = s.substr(0, pos + 1);
			else s.clear();

			if (!s.empty())
				list.push_back(GURL::Filename::UTF8(s.c_str()));

			s.clear();

			if (c < 0) break;
		} else {
			s.push_back(c);
		}
	}
}

int 
main(int argc, const char **argv)
{
  DJVU_LOCALE;
  if (argc < 3) usage();
  GArray<GUTF8String> dargv(0,argc-1);
  for(int i=0;i<argc;++i)
    dargv[i]=GNativeString(argv[i]);
  G_TRY
    {
	  std::vector<GURL> inputlist;
	  std::string outputname;
      cjb2opts opts;
      // Defaults
      opts.forcedpi = 0;
      opts.dpi = 300;
      opts.losslevel = 0;
      opts.verbose = false;
	  opts.pages_per_dict = 10;
      // Parse options
      for (int i=1; i<argc; i++)
        {
          GUTF8String arg = dargv[i];
          if (arg == "-dpi")
		  {
			  if (i + 1 >= argc) usage();
			  char *end;
			  opts.dpi = opts.forcedpi = strtol(dargv[++i], &end, 10);
			  if (*end || opts.dpi < 25 || opts.dpi>1200)
				  usage();
		  }
		  else if (arg == "-dict")
		  {
			  if (i + 1 >= argc) usage();
			  opts.dict = GURL::Filename::UTF8(dargv[++i]);
		  }
		  else if (arg == "-losslevel")
            {
			  if (i + 1 >= argc) usage();
			  char *end;
              opts.losslevel = strtol(dargv[++i], &end, 10);
              if (*end || opts.losslevel<0 || opts.losslevel>200)
                usage();
            }
          else if (arg == "-lossless")
            opts.losslevel = 0;
          else if (arg == "-lossy")
            opts.losslevel = 100;
          else if (arg == "-clean") // almost deprecated
            opts.losslevel = 1;
          else if (arg == "-verbose" || arg == "-v")
            opts.verbose = true;
		  else if (arg == "-pages-per-dict" || arg == "-p")
		  {
			  if (i + 1 >= argc) usage();
			  char *end;
			  opts.pages_per_dict = strtol(dargv[++i], &end, 10);
			  if (*end) usage();
		  } else if (arg == "-output-dict")
		  {
			  if (i + 1 >= argc) usage();
			  opts.output_dict = (const char*)dargv[++i];
		  } else if (arg == "-f") {
			  if (i + 1 >= argc) usage();
			  add_filelist(GURL::Filename::UTF8(dargv[++i]), inputlist);
		  } else if (arg[0] == '-' && arg[1])
			  usage();
		  else if (i == argc - 1) {
			  outputname = (const char*)arg;
		  } else {
			  inputlist.push_back(GURL::Filename::UTF8(arg));
		  }
        }
	  if (inputlist.size() == 0 || outputname.length() == 0)
		  usage();
      // Execute
	  cjb2(inputlist, outputname, opts);
    }
  G_CATCH(ex)
    {
      ex.perror();
      exit(1);
    }
  G_ENDCATCH;
  return 0;
}

