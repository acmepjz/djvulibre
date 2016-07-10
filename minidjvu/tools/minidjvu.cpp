/*
 * minidjvu.cpp - an example of using the library
 */

#include <minidjvu/minidjvu.h>
#include "../src/base/mdjvucfg.h" /* for i18n, HAVE_TIFF */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <locale.h>

/* TODO: remove duplicated code */

#include <string>
#include <vector>
#include <set>

struct FileAndPageNumber {
	std::string fileName;
	int pageNumber;
};

#if HAVE_FREEIMAGE
std::string lastTIFFFileName;
FIMULTIBITMAP *lastTIFF = NULL;
#endif

/* options */
int32 dpi = 300;
int32 pages_per_dict = 10; /* 0 means infinity */
int dpi_specified = 0;
int verbose = 0;
int smooth = 0;
int averaging = 0;
int match = 0;
int Match = 0;
int aggression = 100;
int erosion = 0;
int clean = 0;
int report = 0;
int no_prototypes = 0;
int warnings = 0;
int indirect = 0;
const char* dict_suffix = NULL;
std::vector<FileAndPageNumber> filelist;
int showpage = 0;

/* ========================================================================= */

/* file name template routines (for multipage encoding) {{{ */

static int get_ext_delim_pos(const char *fname)
{
    int pos = strcspn(fname,".");
    int last = 0;
    
    while (last + pos != strlen(fname))
    {
        last += (pos + 1);
        pos = strcspn(fname + last,".");
    }
    return last;
}

static std::string get_page_or_dict_name(const std::set<std::string>& elements2, int cnt, const char *fname, const char *ext)
{
    int extpos;
    std::string page_name;
    
    extpos = get_ext_delim_pos(fname);
	if (extpos > 0)
		page_name = std::string(fname, fname + (extpos - 1));
	else
		page_name = fname;

	if (elements2.find((page_name + ".") + ext) == elements2.end()) {
		return (page_name + ".") + ext;
	}

	for (;; cnt++) {
		char s[16];
		sprintf(s, "#%03d.", cnt);
		if (elements2.find((page_name + s) + ext) == elements2.end()) {
			return (page_name + s) + ext;
		}
	}
}

static void replace_suffix(char *name, const char *suffix)
{
    int len = strlen(name);
    
    name[len-4] = '\0';
    strcat(name, suffix);
}

/* file name template routines (for multipage encoding) }}} */

/* ========================================================================= */

static void show_usage_and_exit(void)           /* {{{ */
{
    const char *what_it_does = _("encode/decode bitonal DjVu files");
    if (strcmp(MDJVU_VERSION, mdjvu_get_version()))
    {
        printf(_("minidjvu - %s\n"), what_it_does);
        printf(_("Warning: program and library version mismatch:\n"));
        printf(_("    program version %s, library version %s.\n\n"), MDJVU_VERSION, mdjvu_get_version());

    }
    else
    {
        printf("minidjvu %s - %s\n", MDJVU_VERSION, what_it_does);

    }
    printf(_("Usage:\n"));
    printf(_("single page encoding/decoding:\n"));
    printf(_("    minidjvu [options] <input file> <output file>\n"));
    printf(_("multiple pages encoding:\n"));
    printf(_("    minidjvu [options] <input file> ... <output file>\n"));
    printf(_("Formats supported:\n"));

#if HAVE_FREEIMAGE
	printf(_("    Any FreeImage-supported format (1-bit only), e.g. BMP, TIF, JPG, PNG, etc.\n"));
#else
    printf(_("    DjVu (single-page bitonal), PBM, Windows BMP"));
    if (mdjvu_have_tiff_support())
        printf(_(", TIFF.\n"));
    else
        printf(_("; TIFF support is OFF.\n"));
#endif

    printf(_("Options:\n"));
    printf(_("    -A, --Averaging:               compute \"average\" representatives\n"));
    printf(_("    -a <n>, --aggression <n>:      set aggression level (default 100)\n"));
    printf(_("    -c, --clean                    remove small black pieces\n"));
    printf(_("    -d <n>, --dpi <n>:             set resolution in dots per inch\n"));
    printf(_("    -e, --erosion                  sacrifice quality to gain in size\n"));
    printf(_("    -f, --filelist <filelist>      input a list of files from a text file\n"));
    printf(_("    -i, --indirect:                generate an indirect multipage document\n"));
    printf(_("    -l, --lossy:                   use all lossy options (-s -c -m -e -A)\n"));
    printf(_("    -L, --Lossy:                   experimental (-s -c -m -M -e -A)\n"));
    printf(_("    -m, --match:                   match and substitute patterns\n"));
    printf(_("    -M, --Match:                   experimental\n"));
    printf(_("    -n, --no-prototypes:           do not search for prototypes\n"));
    printf(_("    -p <n>, --pages-per-dict <n>:  pages per dictionary (default 10)\n"));
    printf(_("    -r, --report:                  report multipage coding progress\n"));
    printf(_("    -R, --Report:                  report progress in the real time\n"));
    printf(_("                                   (good for displaying the progress in a GUI)\n"));
    printf(_("    -s, --smooth:                  remove some badly looking pixels\n"));
    printf(_("    -v, --verbose:                 print messages about everything\n"));
    printf(_("    -X, --Xtension:                file extension for shared dictionary files\n"));
    printf(_("    -w, --warnings:                do not suppress TIFF warnings\n"));
    printf(_("See the man page for detailed description of each option.\n"));
    exit(2);
}                   /* }}} */

static int decide_if_bmp(const char *path)
{
    return mdjvu_ends_with_ignore_case(path, ".bmp");
}

static int decide_if_djvu(const char *path)
{
    return mdjvu_ends_with_ignore_case(path, ".djvu")
        || mdjvu_ends_with_ignore_case(path, ".djv");
}

static int decide_if_tiff(const char *path)
{
#if HAVE_FREEIMAGE
	FREE_IMAGE_FORMAT fif = FreeImage_GetFileType(path, 0);

	if (fif == FIF_TIFF) return 1;
	else return 0;
#else
    return mdjvu_ends_with_ignore_case(path, ".tiff")
        || mdjvu_ends_with_ignore_case(path, ".tif");
#endif
}

/* ========================================================================= */

static mdjvu_image_t load_image(const char *path)
{
    mdjvu_error_t error;
    mdjvu_image_t image;

    if (verbose) printf(_("loading a DjVu page from `%s'\n"), path);
    image = mdjvu_load_djvu_page(path, &error);
    if (!image)
    {
        fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
        exit(1);
    }
    if (verbose)
    {
        printf(_("loaded; the page has %d bitmaps and %d blits\n"),
               mdjvu_image_get_bitmap_count(image),
               mdjvu_image_get_blit_count(image));
    }
    return image;
}

static mdjvu_matcher_options_t get_matcher_options(void)
{
    mdjvu_matcher_options_t m_options = NULL;
    if (match || Match)
    {
        m_options = mdjvu_matcher_options_create();
        mdjvu_use_matcher_method(m_options, MDJVU_MATCHER_PITH_2);
        if (Match)
            mdjvu_use_matcher_method(m_options, MDJVU_MATCHER_RAMPAGE);
        mdjvu_set_aggression(m_options, aggression);
    }
    return m_options;
}

static void sort_and_save_image(mdjvu_image_t image, const char *path)
{
    mdjvu_error_t error;

    mdjvu_compression_options_t options = mdjvu_compression_options_create();
    mdjvu_set_matcher_options(options, get_matcher_options());

    mdjvu_set_clean(options, clean);
    mdjvu_set_verbose(options, verbose);
    mdjvu_set_no_prototypes(options, no_prototypes);
    mdjvu_set_averaging(options, averaging);
    mdjvu_compress_image(image, options);
    mdjvu_compression_options_destroy(options);

    if (verbose) printf(_("encoding to `%s'\n"), path);

    if (!mdjvu_save_djvu_page(image, path, NULL, &error, erosion))
    {
        fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
        exit(1);
    }
}

static mdjvu_bitmap_t load_bitmap(const char *path, int tiff_idx)
{
    mdjvu_error_t error;
    mdjvu_bitmap_t bitmap;

#if HAVE_FREEIMAGE
	if (lastTIFFFileName == path && lastTIFF != NULL)
	{
		FIBITMAP *dib = FreeImage_LockPage(lastTIFF, tiff_idx);
		if (dib)
		{
			bitmap = mdjvu_file_load_fibitmap(dib, &error);
			FreeImage_UnlockPage(lastTIFF, dib, FALSE);
		}
		if (!bitmap)
		{
			fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
			exit(1);
		}
	}
	else
#endif
    if (decide_if_bmp(path))
    {
        if (verbose) printf(_("loading from Windows BMP file `%s'\n"), path);
        bitmap = mdjvu_load_bmp(path, &error);
    }
    else if (decide_if_tiff(path))
    {
#if HAVE_FREEIMAGE
		if (lastTIFF) {
			FreeImage_CloseMultiBitmap(lastTIFF, 0);
			lastTIFF = NULL;
		}
		lastTIFF = FreeImage_OpenMultiBitmap(FIF_TIFF, path, FALSE, TRUE, TRUE, 0);
		FIBITMAP *dib = FreeImage_LockPage(lastTIFF, tiff_idx);
		if (dib)
		{
			bitmap = mdjvu_file_load_fibitmap(dib, &error);
			FreeImage_UnlockPage(lastTIFF, dib, FALSE);
		}
		if (!bitmap)
		{
			fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
			exit(1);
		}
#elif HAVE_TIFF
        if (verbose) printf(_("loading from TIFF file `%s'\n"), path);
        if (!warnings)
            mdjvu_disable_tiff_warnings();
        if (dpi_specified)
            bitmap = mdjvu_load_tiff(path, NULL, &error, tiff_idx);
        else
            bitmap = mdjvu_load_tiff(path, &dpi, &error, tiff_idx);
        if (verbose) printf(_("resolution is %d dpi\n"), dpi);
#endif
    }
    else if (decide_if_djvu(path))
    {
        mdjvu_image_t image = load_image(path);
        bitmap = mdjvu_render(image);
        mdjvu_image_destroy(image);
        if (verbose)
        {
            printf(_("bitmap %d x %d rendered\n"),
                   mdjvu_bitmap_get_width(bitmap),
                   mdjvu_bitmap_get_height(bitmap));
        }
    }
    else
    {
#if HAVE_FREEIMAGE
		if (verbose) printf(_("loading from graphical file `%s'\n"), path);
		bitmap = mdjvu_load_fibitmap(path, &error);
#else
        if (verbose) printf(_("loading from PBM file `%s'\n"), path);
        bitmap = mdjvu_load_pbm(path, &error);
#endif
    }

    if (!bitmap)
    {
        fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
        exit(1);
    }

    if (smooth)
    {
        if (verbose) printf(_("smoothing the bitmap\n"));
        mdjvu_smooth(bitmap);
    }

    return bitmap;
}

static void save_bitmap(mdjvu_bitmap_t bitmap, const char *path)
{
    mdjvu_error_t error;
    int result;

    if (decide_if_bmp(path))
    {
        if (verbose) printf(_("saving to Windows BMP file `%s'\n"), path);
        result = mdjvu_save_bmp(bitmap, path, dpi, &error);
    }
    else if (decide_if_tiff(path))
    {
        if (verbose) printf(_("saving to TIFF file `%s'\n"), path);
        if (!warnings)
            mdjvu_disable_tiff_warnings();
        result = mdjvu_save_tiff(bitmap, path, &error);
    }
    else
    {
        if (verbose) printf(_("saving to PBM file `%s'\n"), path);
        result = mdjvu_save_pbm(bitmap, path, &error);
    }

    if (!result)
    {
        fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
        exit(1);
    }
}

/* ========================================================================= */

static void decode(const char* infile, const char* outfile)
{
    mdjvu_image_t image;    /* a sequence of blits (what is stored in DjVu) */
    mdjvu_bitmap_t bitmap;  /* the result                                   */

    if (verbose) printf(_("\nDECODING\n"));
    if (verbose) printf(_("________\n\n"));

	image = load_image(infile);
    bitmap = mdjvu_render(image);
    mdjvu_image_destroy(image);

    if (verbose)
    {
        printf(_("bitmap %d x %d rendered\n"),
               mdjvu_bitmap_get_width(bitmap),
               mdjvu_bitmap_get_height(bitmap));
    }

    if (smooth)
    {
        if (verbose) printf(_("smoothing the bitmap\n"));
        mdjvu_smooth(bitmap);
    }

	save_bitmap(bitmap, outfile);
    mdjvu_bitmap_destroy(bitmap);
}


static mdjvu_image_t split_and_destroy(mdjvu_bitmap_t bitmap)
{
    mdjvu_image_t image;
    if (verbose) printf(_("splitting the bitmap into pieces\n"));
    image = mdjvu_split(bitmap, dpi, /* options:*/ NULL);
    mdjvu_bitmap_destroy(bitmap);
    if (verbose)
    {
        printf(_("the splitted image has %d pieces\n"),
                mdjvu_image_get_blit_count(image));
    }
    if (clean)
    {
        if (verbose) printf(_("cleaning\n"));
        mdjvu_clean(image);
        if (verbose)
        {
            printf(_("the cleaned image has %d pieces\n"),
                    mdjvu_image_get_blit_count(image));
        }
    }
    return image;
}


static void encode(const FileAndPageNumber& infile, const char* outfile)
{
    mdjvu_bitmap_t bitmap;
    mdjvu_image_t image;

    if (verbose) printf(_("\nENCODING\n"));
    if (verbose) printf(_("________\n\n"));

	bitmap = load_bitmap(infile.fileName.c_str(), infile.pageNumber);

    image = split_and_destroy(bitmap);
	sort_and_save_image(image, outfile);
    mdjvu_image_destroy(image);
}


/* Filtering is nondjvu->nondjvu job. */
static void filter(const FileAndPageNumber& infile, const char* outfile)
{
    mdjvu_bitmap_t bitmap;

    if (verbose) printf(_("\nFILTERING\n"));
    if (verbose) printf(_("_________\n\n"));

	bitmap = load_bitmap(infile.fileName.c_str(), infile.pageNumber);
	save_bitmap(bitmap, outfile);
    mdjvu_bitmap_destroy(bitmap);
}


static const char *strip(const char *str, char sep)
{
    const char *t = strrchr(str, sep);
    if (t)
        return t + 1;
    else
        return str;
}

/* return path without a directory name */ 
static const char *strip_dir(const char *path)
{
    return strip(strip(path, '\\'), '/');
}


static void multipage_encode(const std::vector<FileAndPageNumber> &pages, const char *outname)
{
    mdjvu_image_t *images;
    mdjvu_image_t dict;
	const int n = pages.size();
    int i, el = 0;
    const int ndicts = (pages_per_dict <= 0)? 1 :
                                        (n % pages_per_dict > 0) ?  (int) (n/pages_per_dict) + 1:
                                                                    (int) (n/pages_per_dict);
	std::string dict_name, path;
	std::vector<std::string> elements(n + ndicts);
	std::set<std::string> elements2;
	std::vector<int> sizes(n + ndicts, 0);
    mdjvu_compression_options_t options;
    mdjvu_bitmap_t bitmap;
    mdjvu_error_t error;
    int32 pages_compressed;
    FILE *f, *tf=NULL;

    match = 1;

    if (!decide_if_djvu(outname))
    {
        fprintf(stderr, _("when encoding many pages, output file must be DjVu\n"));
        exit(1);
    }
    if (!indirect)
    {
        tf = tmpfile();
        if (!tf)
        {
            fprintf(stderr, _("Could not create a temporary file\n"));
            exit(1);
        }
    }

    if (verbose) printf(_("\nMULTIPAGE ENCODING\n"));
    if (verbose) printf(_("__________________\n\n"));
    if (verbose) printf(_("%d pages total\n"), n);

    options = mdjvu_compression_options_create();
    mdjvu_set_matcher_options(options, get_matcher_options());

    mdjvu_set_clean(options, clean);
    mdjvu_set_verbose(options, verbose);
    mdjvu_set_no_prototypes(options, no_prototypes);
    mdjvu_set_report(options, report);
	mdjvu_set_showpage(options, showpage);
    mdjvu_set_averaging(options, averaging);
    mdjvu_set_report_total_pages(options, n);

    /* compressing */
    if (pages_per_dict <= 0) pages_per_dict = n;
    if (pages_per_dict > n) pages_per_dict = n;
    images = MDJVU_MALLOCV(mdjvu_image_t, pages_per_dict);
    pages_compressed = 0;

    while (n - pages_compressed)
    {
        int32 pages_to_compress = n - pages_compressed;
        if (pages_to_compress > pages_per_dict)
            pages_to_compress = pages_per_dict;

        mdjvu_set_report_start_page(options, pages_compressed + 1);

        for (i = 0; i < pages_to_compress; i++)
        {
			bitmap = load_bitmap(pages[pages_compressed + i].fileName.c_str(), pages[pages_compressed + i].pageNumber);
            images[i] = split_and_destroy(bitmap);
            if (report)
                printf(_("Loading: %d of %d completed\n"), pages_compressed + i + 1, n);
        }

        dict = mdjvu_compress_multipage(pages_to_compress, images, options);

		dict_name = get_page_or_dict_name(elements2, pages[pages_compressed].pageNumber, strip_dir(pages[pages_compressed].fileName.c_str()), dict_suffix);
        
        if (!indirect)
            sizes[el] = mdjvu_file_save_djvu_dictionary(dict, (mdjvu_file_t) tf, 0, &error, erosion);
        else
            sizes[el] = mdjvu_save_djvu_dictionary(dict, dict_name.c_str(), &error, erosion);
        
        if (!sizes[el])
        {
            fprintf(stderr, "%s: %s\n", dict_name.c_str(), mdjvu_get_error_message(error));
            exit(1);
        }
        elements[el++] = dict_name;
		elements2.insert(dict_name);

        for (i = 0; i < pages_to_compress; i++)
        {
			path = get_page_or_dict_name(elements2, pages[pages_compressed + i].pageNumber,
				strip_dir(pages[pages_compressed + i].fileName.c_str()), "djvu");

            if (verbose)
                printf(_("saving page #%d into %s using dictionary %s\n"), pages_compressed + i + 1, path.c_str(), dict_name.c_str());
            
            if (!indirect)
                sizes[el] = mdjvu_file_save_djvu_page(images[i], (mdjvu_file_t) tf, strip_dir(dict_name.c_str()), 0, &error, erosion);
            else
                sizes[el] = mdjvu_save_djvu_page(images[i], path.c_str(), strip_dir(dict_name.c_str()), &error, erosion);
            if (!sizes[el])
            {
                fprintf(stderr, "%s: %s\n", path.c_str(), mdjvu_get_error_message(error));
                exit(1);
            }
            elements[el++] = path;
			elements2.insert(path);
            mdjvu_image_destroy(images[i]);
            if (report)
                printf(_("Saving: %d of %d completed\n"), pages_compressed + i + 1, n);
        }
        mdjvu_image_destroy(dict);
        pages_compressed += pages_to_compress;
    }

	std::vector<const char*> elements0;
	for (i = 0; i < el; i++) {
		elements0.push_back(elements[i].c_str());
	}

	if (!indirect)
	{
		f = fopen(outname, "wb");
		if (!f)
		{
			fprintf(stderr, "%s: %s\n", outname, (const char *)mdjvu_get_error(mdjvu_error_fopen_write));
			exit(1);
		}
		mdjvu_file_save_djvu_dir(const_cast<char**>(&(elements0[0])), &(sizes[0]), el, (mdjvu_file_t)f, (mdjvu_file_t)tf, &error);
		fclose(tf);
		fclose(f);
	} else
		mdjvu_save_djvu_dir(const_cast<char**>(&(elements0[0])), &(sizes[0]), el, outname, &error);

    /* destroying */
    mdjvu_compression_options_destroy(options);

    MDJVU_FREEV(images);
}

/* same_option(foo, "opt") returns 1 in three cases:
 *
 *      foo is "o" (first letter of opt)
 *      foo is "opt"
 *      foo is "-opt"
 */
static int same_option(const char *option, const char *s)
{
    if (option[0] == s[0] && !option[1]) return 1;
    if (!strcmp(option, s)) return 1;
    if (option[0] == '-' && !strcmp(option + 1, s)) return 1;
    return 0;
}

static void add_file(const char* filename) {
	int pagecount = 1;

	if (decide_if_tiff(filename)) {
		pagecount = mdjvu_get_tiff_page_count(filename);
	}

	FileAndPageNumber page = { filename, 0 };
	for (int i = 0; i < pagecount; i++) {
		page.pageNumber = i;
		filelist.push_back(page);
	}
}

static void add_filelist(const char* filename) {
	FILE *f = fopen(filename, "rb");
	if (f == NULL) {
		fprintf(stderr, "Can't load file list '%s'\n", filename);
		exit(2);
	}

	int c;
	std::string s;

	for (;;) {
		c = fgetc(f);

		if (c == '\r') continue;
		else if (c < 0 || c == '\n') {
			size_t pos = s.find_first_not_of(" \t");
			if (pos != s.npos) s = s.substr(pos);
			else s.clear();
			pos = s.find_last_not_of(" \t");
			if (pos != s.npos) s = s.substr(0, pos + 1);
			else s.clear();

			add_file(s.c_str());

			s.clear();

			if (c < 0) break;
		} else {
			s.push_back(c);
		}
	}

	fclose(f);
}

static int process_options(int argc, char **argv)
{
    int i;
    for (i = 1; i < argc; i++)
    {
		if (argv[i][0] != '-') {
			add_file(argv[i]);
			continue;
		}

		char *option = argv[i] + 1;
        if (same_option(option, "verbose"))
            verbose = 1;
        else if (same_option(option, "smooth"))
            smooth = 1;
        else if (same_option(option, "match"))
            match = 1;
        else if (same_option(option, "Match"))
            Match = 1;
        else if (same_option(option, "no-prototypes"))
            no_prototypes = 1;
        else if (same_option(option, "erosion"))
            erosion = 1;
        else if (same_option(option, "clean"))
            clean = 1;
        else if (same_option(option, "warnings"))
            warnings = 1;
		else if (same_option(option, "report"))
			report = 1;
		else if (same_option(option, "Report"))
			showpage = 1;
        else if (same_option(option, "Averaging"))
            averaging = 1;
        else if (same_option(option, "lossy"))
        {
            smooth = 1;
            match = 1;
            erosion = 1;
            clean = 1;
            averaging = 1;
        }
        else if (same_option(option, "Lossy"))
        {
            smooth = 1;
            Match = match = 1;
            erosion = 1;
            clean = 1;
            averaging = 1;
        }
        else if (same_option(option, "pages-per-dict"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            pages_per_dict = atoi(argv[i]);
            if (pages_per_dict < 0)
            {
                fprintf(stderr, _("bad --pages-per-dict value\n"));
                exit(2);
            }
        }
        else if (same_option(option, "dpi"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            dpi = atoi(argv[i]);
            dpi_specified = 1;
            if (dpi < 20 || dpi > 2000)
            {
                fprintf(stderr, _("bad resolution\n"));
                exit(2);
            }
        }
		else if (same_option(option, "filelist"))
		{
			i++;
			if (i == argc) show_usage_and_exit();
			add_filelist(argv[i]);
		}
		else if (same_option(option, "aggression"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            aggression = atoi(argv[i]);

            match = 1;
        }
        else if (same_option(option, "Xtension"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            dict_suffix = argv[i];
        }
        else if (same_option(option, "indirect"))
            indirect = 1;
        else
        {
            fprintf(stderr, _("unknown option: %s\n"), argv[i]);
            exit(2);
        }
    }
    return i;
}

int main(int argc, char **argv)
{
    setlocale(LC_ALL, "");
#ifdef HAVE_GETTEXT
    bindtextdomain("minidjvu", LOCALEDIR);
    textdomain("minidjvu");
#endif


	process_options(argc - 1, argv);
    if ( dict_suffix == NULL ) dict_suffix = "iff";

    if (argc < 3 || filelist.empty())
        show_usage_and_exit();

	if (filelist.size() >= 2)
    {
		multipage_encode(filelist, argv[argc - 1]);
    }
    else if (decide_if_djvu(argv[argc - 1]))
    {
		encode(filelist[0], argv[argc - 1]);
    }
    else
    {
		if (decide_if_djvu(filelist[0].fileName.c_str()))
			decode(filelist[0].fileName.c_str(), argv[argc - 1]);
        else
			filter(filelist[0], argv[argc - 1]);
    }

#if HAVE_FREEIMAGE
	if (lastTIFF) {
		FreeImage_CloseMultiBitmap(lastTIFF, 0);
		lastTIFF = NULL;
	}
#endif

    if (verbose) printf("\n");
    #ifndef NDEBUG 
        if (alive_bitmap_counter)
           printf(_("alive_bitmap_counter = %d\n"), alive_bitmap_counter);
    #endif

    return 0;
}
