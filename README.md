# rts2nii

This Python project is designed to convert DICOM RT Structure files (saved in form of long sequence of sequences per ROI) containing 2D segmentation masks into NIfTI format.

## Dependencies

- pydicom
- nibabel
- numpy
- matplotlib
- skimage
- os
more details in [requirements.txt](https://github.com/xkuubix/rts2nii/blob/main/requirements.txt)

## DICOM Tags

The script uses several [DICOM tags](https://dicom.innolitics.com/ciods) for linking, identification, image properties, referencing, sequences, and ROI and contour sequences. These tags are defined as constants at the top of the script.

## Usage

To use this script, one need to have a DICOM RT Structure file with connected DICOM images to convert masks into NIfTI format. Script iterates over every directory in given path. Masks are saved in folder containing original mask in DICOM_RT-struct. Original file name is joined with ROI_NAME tags from RT-struct that are present on saved mask.

## Future Work

The project is still under development. Future updates will include more file formats, handling overlaying ROIs, error handling, CLI interface.

## Contributing

Contributions are welcome. Please open an issue to discuss your idea or submit a Pull Request with your changes.

## License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/xkuubix/rts2nii/blob/main/LICENSE) file for details.
