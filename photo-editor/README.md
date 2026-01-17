# Photo Editor Web Application

A simple local web application for viewing, cropping, and rotating JPEG photos one at a time.

## Features

- View photos from any directory on your computer
- Navigate through photos one at a time (Previous/Next)
- Rotate photos left or right (90-degree increments)
- Crop photos with an intuitive drag-and-drop interface
- Save changes directly back to the original files
- Keyboard shortcuts for quick navigation
- Clean, modern interface

## Installation

1. Make sure you have Node.js installed (v12 or higher)
2. Navigate to this directory
3. Dependencies are already installed

## Usage

### Starting the Application

1. Open a terminal in this directory
2. Start the server:
   ```bash
   node server.js
   ```

   Or specify a photo directory:
   ```bash
   PHOTO_DIR=/path/to/your/photos node server.js
   ```

3. Open your browser and go to: `http://localhost:3000`

### Using the Interface

1. **Set Photo Directory**
   - Enter the full path to your photo directory in the text field at the top
   - Click "Set Directory" to change the directory
   - Click "Load Photos" to reload images from the current directory

2. **Navigate Photos**
   - Use "Previous" and "Next" buttons to move through photos
   - Or use left/right arrow keys on your keyboard
   - Click on any image name in the sidebar to jump to it

3. **Rotate Photos**
   - Click "Rotate Left" to rotate 90° counter-clockwise
   - Click "Rotate Right" to rotate 90° clockwise

4. **Crop Photos**
   - Click "Enable Crop" to enter crop mode
   - Drag the corners or edges of the crop box to adjust
   - Move the entire crop box by dragging inside it
   - Click "Reset Crop" to start over
   - Click "Disable Crop" to exit crop mode

5. **Save Changes**
   - Click "Save Changes" (or press Ctrl/Cmd+S) to save
   - The original file will be overwritten with your edits
   - A success message will appear when saved

6. **Reset**
   - Click "Reset Image" to discard all unsaved changes

### Keyboard Shortcuts

- **Left Arrow**: Previous image
- **Right Arrow**: Next image
- **Ctrl+S** / **Cmd+S**: Save current image

## Technical Details

### Supported Formats
- JPEG (.jpg, .jpeg)

### Technologies Used
- **Backend**: Node.js, Express, Sharp (image processing)
- **Frontend**: HTML, CSS, JavaScript, Cropper.js
- **Port**: 3000 (default)

### File Structure
```
photo-editor/
├── server.js           # Backend server
├── package.json        # Dependencies
├── sample-photos/      # Default photo directory (empty)
├── public/
│   ├── index.html     # Main interface
│   ├── style.css      # Styling
│   └── app.js         # Frontend logic
└── README.md          # This file
```

## Important Notes

⚠️ **WARNING**: This application saves changes directly to your original files. Make sure you have backups of important photos before editing!

- Changes are permanent once saved
- No undo functionality (except resetting before saving)
- The application modifies files in place

## Troubleshooting

**Port already in use**
- If port 3000 is already taken, modify the `PORT` variable in `server.js`

**Images not loading**
- Make sure the directory path is correct and accessible
- Check that the directory contains JPEG files
- Verify file permissions allow reading/writing

**Can't save changes**
- Ensure you have write permissions for the photo directory
- Check that the files aren't read-only

## Example Usage

```bash
# Start with default directory (sample-photos)
node server.js

# Start with specific directory
PHOTO_DIR=/home/user/Pictures node server.js

# Or set the directory through the web interface
```

## Security Notes

- This application is intended for LOCAL USE ONLY
- Do not expose it to the internet
- Only access photos from trusted directories
- The app includes basic path traversal protection

## License

Free to use and modify as needed.
