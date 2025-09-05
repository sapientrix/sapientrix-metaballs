import React from 'react'
import { createRoot } from 'react-dom/client'
import MetaballsGrid from './metaballs'
import './index.css'

const rootEl = document.getElementById('root')!
const root = createRoot(rootEl)
root.render(<MetaballsGrid />)
