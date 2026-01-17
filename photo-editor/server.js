const express = require('express');
const multer = require('multer');
const sharp = require('sharp');
const path = require('path');
const fs = require('fs').promises;
const cors = require('cors');

const app = express();
const PORT = 3000;

app.use(cors());
app.use(express.json({ limit: '50mb' }));
app.use(express.static('public'));

// Configuration - change this to your photo directory
let PHOTO_DIR = process.env.PHOTO_DIR || path.join(__dirname, 'sample-photos');

// Get list of images in directory
app.get('/api/images', async (req, res) => {
  try {
    const files = await fs.readdir(PHOTO_DIR);
    const imageFiles = files.filter(file => {
      const ext = path.extname(file).toLowerCase();
      return ['.jpg', '.jpeg'].includes(ext);
    });
    res.json({ images: imageFiles, directory: PHOTO_DIR });
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// Get a specific image
app.get('/api/image/:filename', async (req, res) => {
  try {
    const filename = req.params.filename;
    const filepath = path.join(PHOTO_DIR, filename);

    // Check if file exists and is within PHOTO_DIR (security check)
    const realPath = await fs.realpath(filepath);
    const realPhotoDir = await fs.realpath(PHOTO_DIR);

    if (!realPath.startsWith(realPhotoDir)) {
      return res.status(403).json({ error: 'Access denied' });
    }

    res.sendFile(realPath);
  } catch (error) {
    res.status(404).json({ error: 'Image not found' });
  }
});

// Save edited image
app.post('/api/save', async (req, res) => {
  try {
    const { filename, imageData, rotation, cropData } = req.body;

    if (!filename || !imageData) {
      return res.status(400).json({ error: 'Missing filename or image data' });
    }

    const filepath = path.join(PHOTO_DIR, filename);

    // Remove data URL prefix
    const base64Data = imageData.replace(/^data:image\/\w+;base64,/, '');
    const buffer = Buffer.from(base64Data, 'base64');

    // Process image with sharp
    let image = sharp(buffer);

    // Apply rotation if specified
    if (rotation && rotation !== 0) {
      image = image.rotate(rotation);
    }

    // Save the image
    await image.toFile(filepath);

    res.json({ success: true, message: 'Image saved successfully' });
  } catch (error) {
    console.error('Error saving image:', error);
    res.status(500).json({ error: error.message });
  }
});

// Update photo directory
app.post('/api/set-directory', async (req, res) => {
  try {
    const { directory } = req.body;

    if (!directory) {
      return res.status(400).json({ error: 'Directory path required' });
    }

    // Check if directory exists
    const stats = await fs.stat(directory);
    if (!stats.isDirectory()) {
      return res.status(400).json({ error: 'Path is not a directory' });
    }

    PHOTO_DIR = directory;
    res.json({ success: true, directory: PHOTO_DIR });
  } catch (error) {
    res.status(400).json({ error: 'Directory not found or inaccessible' });
  }
});

app.listen(PORT, () => {
  console.log(`Photo Editor running at http://localhost:${PORT}`);
  console.log(`Photo directory: ${PHOTO_DIR}`);
  console.log(`\nTo change photo directory, set PHOTO_DIR environment variable:`);
  console.log(`  PHOTO_DIR=/path/to/photos node server.js`);
});
