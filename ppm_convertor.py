
from PIL import Image

# Open the PPM file
with Image.open('ukr.ppm') as img:
    # Save it as a JPEG file
    img.save('ukr.jpg', 'JPEG')

print("Conversion completed successfully.")
