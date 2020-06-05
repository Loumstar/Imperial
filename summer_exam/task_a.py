from matplotlib import pyplot as plot

banksy_filename = "Banksy.jpg"

# Read and show image
banksy_image = plot.imread(banksy_filename).copy()
plot.imshow(banksy_image)
plot.show()

def create_blue_red_triangles(image):
    """
    Method to filter an image such that the upper triangle only has blue colour and the
    lower triangle only has red colour.

    Returns a copy of the input image with the changes.
    """
    # Determine image dimensions
    m, n = len(image), len(image[0])
    
    # Copy image so as to not affect the original
    triangle_image = image.copy()

    for j, row in enumerate(image):
        for i, _ in enumerate(row):
            # If the coordinate i,j exists in the upper triangle
            if abs(i - (n / 2)) / n * m < (m / 2) - j:
                triangle_image[j][i][0] = 0 # set red to 0
                triangle_image[j][i][1] = 0 # set green to 0
            # If the coordinate i,j exists in the lower triangle
            elif abs(i - (n / 2)) / n * m < j - (m / 2):
                triangle_image[j][i][1] = 0 # set green to 0
                triangle_image[j][i][2] = 0 # set blue to 0

    return triangle_image

# Save the new image
plot.imsave('Tactics.jpg', create_blue_red_triangles(banksy_image))