#include <serror.h>
#include <spng.h>
#include <png.h>

// Newer versions of libpng do not define these
#ifndef png_infopp_NULL
#define png_infopp_NULL (png_infopp)NULL
#define int_p_NULL (int*)NULL
#endif

namespace skn
{
	// Write a 24-bit rgb png in (x,y) ordering
	void write_png(const char * filename, const Image & img)
	{
		FILE * file = fopen(filename,"w");
		png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
		if(!png_ptr) serror("Could not allocate png write structure!");
		png_infop info_ptr = png_create_info_struct(png_ptr);
		if(!info_ptr) serror("Could not allocate png info structure!");
		png_bytep * row_pointers = 0;
		if(setjmp(png_jmpbuf(png_ptr)))
		{
			png_destroy_write_struct(&png_ptr, &info_ptr);
			if(row_pointers) {
				for(int i = 0; i < img.size(); i++)
					png_free(png_ptr, row_pointers[i]);
				png_free(png_ptr, row_pointers);
			}
			fclose(file);
			serror("Internal png error!");
		}
		png_init_io(png_ptr, file);

		png_set_IHDR(png_ptr, info_ptr, img.size(0), img.size(1), 8,
				PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
				PNG_FILTER_TYPE_DEFAULT);

		row_pointers = (png_bytep*) png_malloc(png_ptr, img.size(1)*sizeof(png_bytep));
		int pixel_size = 3;
		for(int i = 0; i < img.size(1); i++)
			row_pointers[i] = (png_bytep) png_malloc(png_ptr, img.size(0)*pixel_size);

		// Whew! After all this messy setup, we are finally ready to fill
		// the image structure with data.
		for(int i = 0; i < img.size(0); i++)
			for(int j = 0; j < img.size(1); j++)
				for(int k = 0; k < pixel_size; k++)
					row_pointers[j][i*pixel_size+k] = img(i,j)(k);

		// And write
		png_set_rows(png_ptr, info_ptr, row_pointers);
		png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

		for(int i = 0; i < img.size(1); i++)
			png_free(png_ptr, row_pointers[i]);
		png_free(png_ptr, row_pointers);
		png_destroy_write_struct(&png_ptr, &info_ptr);
		fclose(file);
	}

	Image read_png(const char * filename)
	{
		FILE * fp = fopen(filename, "r");
		png_uint_32 width, height;
		int bit_depth, color_type, interlace_type;

		png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

		if (png_ptr == NULL)
		{
			fclose(fp);
			serror("Error creating png read struct!");
		}

		/* Allocate/initialize the memory for image information.  REQUIRED. */
		png_infop info_ptr = png_create_info_struct(png_ptr);
		if (info_ptr == NULL)
		{
			fclose(fp);
			png_destroy_read_struct(&png_ptr, png_infopp_NULL, png_infopp_NULL);
			serror("Error creating png info struct!");
		}

		/* Set error handling if you are using the setjmp/longjmp method (this is
		 * the normal method of doing things with libpng).  REQUIRED unless you
		 * set up your own error handlers in the png_create_read_struct() earlier.
		 */

		if (setjmp(png_jmpbuf(png_ptr)))
		{
			/* Free all of the memory associated with the png_ptr and info_ptr */
			png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
			fclose(fp);
			/* If we get here, we had a problem reading the file */
			serror("Error reading png file!");
		}

		/* Set up the input control if you are using standard C streams */
		png_init_io(png_ptr, fp);
		png_read_info(png_ptr, info_ptr);
		png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
			 &interlace_type, NULL, NULL);

		if(bit_depth == 16) png_set_strip_16(png_ptr);
		else if(bit_depth < 8) png_set_packing(png_ptr);
		if(color_type == PNG_COLOR_TYPE_PALETTE)png_set_expand(png_ptr);
		bool gray = color_type == PNG_COLOR_TYPE_GRAY;

		png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

		if (png_get_color_type(png_ptr, info_ptr) != PNG_COLOR_TYPE_RGB)
		{
			fclose(fp);
			png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
			serror("Image must be 8-bit RBG");
		}

		png_bytep * row_pointers = png_get_rows(png_ptr, info_ptr);
		/* At this point you have read the entire image */
		/* I need the image to be in 24-bit RGB format. It seems
		 * libpng can't do that transformation for me, which is very
		 * irritating. So I will just assme that it is 24-bit RGB anyway */
		Image img(png_get_image_width(png_ptr, info_ptr),
			png_get_image_height(png_ptr, info_ptr));
		int pixel_size = 3;
		for(int i = 0; i < img.size(0); i++)
			for(int j = 0; j < img.size(1); j++)
				for(int k = 0; k < pixel_size; k++)
					img(i,j)(k) = row_pointers[j][i*pixel_size+ (gray ? 1 : k)];
		fclose(fp);
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		return img;
	}
}
