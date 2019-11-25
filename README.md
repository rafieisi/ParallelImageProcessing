# ParallelImageProcessing
# Image Processing
Image processing techniques are exciting and useful in many situations: to obtain various artistic effects, sharpen blurry photos, perform edge detection, etc.

For this assignment you will be working with monochrome (greyscale) images. Monochrome images are represented as a two dimensional array of pixels. Each pixel is stored as a byte and encodes a greyscale colour as a value between 0 and 255.

A discrete Laplace operator is used to compute the second derivatives of an image, which can emphasize edges within an image. This is useful in image processing for performing edge detection and various other related applications. The discrete Laplacian filter is a 3 x 3 array, which typically contains a high negative value at the centre, surrounded by small positive values. Some variations include opposite signs, to achieve a similar effect.

A Laplacian filter can be applied to an existing image, by considering each pixel and its surrounding 8 neighbours and multiplying the 9 pixel values with the corresponding values in the Laplacian filter. The sum of the pairwise multiplications is used as the new value of that pixel.

Pixels on the edges and corners of the image are dealt with a bit differently, since they do not have all 8 neighbours. Only the valid neighbours and the corresponding filter weights are factored into computing the new value of such a "marginal" pixel.

Intuitively, the Laplacian filter emphasizes sudden changes in pixel values. The weights in the Laplacian filter end up setting a darker tone for areas with low pixel changes, and contrasting white tones for sharp changes which represent edges.

Once is pixel value is updated, you might notice that the new pixel values may end up outside the [0,255] range. When processing an image, the new pixel values must be brought back to the original range by a process called normalization. This involves finding the minimum (Min) and maximum (Max), from the new (out-of-bounds) pixel values, then using these to normalize the pixel values back into the [0,255] range. For example, assume that the new pixel values are in the range [-100, 190]. In this case, you must add 100 to all pixels to bring all pixels in the range [0, 290]. Then, you must scale all pixels multiplying with 255/290, to bring all pixels in the [0,255] range.

[Program]
For this assignment, this program has to apply Laplacian filter to each pixel in different methedologies:
    1)Sequential Method
    2)Horizontal Sharding
    3)Vertical Sharding
    4)Work Pool Implementation
