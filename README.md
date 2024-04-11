# rts2nii

This Python project is designed to convert DICOM RT Structure files (saved in form of long sequence of sequences per ROI) containing 2D segmentation masks into NIfTI format.

## Dependencies

- pydicom
- nibabel
- numpy
- matplotlib
- skimage
- os

## DICOM Tags

The script uses several DICOM tags for linking, identification, image properties, referencing, sequences, and ROI and contour sequences. These tags are defined as constants at the top of the script.

## Logging

The script uses Python's built-in logging module to log events that happen while the program runs. The logger is configured to log at the INFO level.

## Usage

To use this script, one need to have a DICOM RT Structure file with connected DICOM images to convert masks into NIfTI format.

## Future Work

The project is still under development. Future updates will include more file formats, handling overlaying ROIs, error handling.

## Contributing

Contributions are welcome. Please open an issue to discuss your idea or submit a Pull Request with your changes.

## License

This project is licensed under the MIT License. See the LICENSE file for details.